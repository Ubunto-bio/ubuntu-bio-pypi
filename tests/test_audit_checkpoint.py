"""Tests for external audit checkpoint RS256 signature and anchor continuity checks."""

from __future__ import annotations

import base64
import hashlib
import json
import subprocess
import time
from pathlib import Path
from typing import Any

import pytest
from typer.testing import CliRunner

from ubuntu_bio.audit import (
    fetch_and_verify_checkpoint,
    load_anchor_ledger,
    verify_checkpoint_continuity,
    verify_rs256_checkpoint_jwt,
)
from ubuntu_bio.cli import app
from ubuntu_bio.exceptions import UbuntuBioError


def _b64url(data: bytes) -> str:
    return base64.urlsafe_b64encode(data).rstrip(b"=").decode("ascii")


def _generate_rsa_key_via_openssl(key_path: Path) -> tuple[int, int]:
    """Generate a 2048-bit RSA keypair with system openssl and return (n, e)."""
    subprocess.run(
        ["openssl", "genrsa", "-out", str(key_path), "2048"],
        check=True,
        capture_output=True,
    )
    mod_out = subprocess.run(
        ["openssl", "rsa", "-in", str(key_path), "-modulus", "-noout"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    n = int(mod_out.split("=", 1)[1], 16)
    return n, 65537


def _sign_test_checkpoint_jwt(
    *,
    key_dir: Path,
    repo: str,
    head_hash: str,
    entries_checked: int,
    kid: str = "test-kid-1",
    key_path: Path | None = None,
) -> tuple[str, dict[str, Any], Path]:
    """Sign a real RS256 JWT using system openssl (zero extra Python crypto dependencies)."""
    pem_file = key_path or (key_dir / "test_audit.pem")
    if not pem_file.exists():
        _generate_rsa_key_via_openssl(pem_file)
    mod_out = subprocess.run(
        ["openssl", "rsa", "-in", str(pem_file), "-modulus", "-noout"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    n = int(mod_out.split("=", 1)[1], 16)
    e = 65537
    n_bytes = n.to_bytes((n.bit_length() + 7) // 8, "big")
    e_bytes = e.to_bytes((e.bit_length() + 7) // 8, "big")
    jwks: dict[str, Any] = {
        "keys": [
            {
                "kty": "RSA",
                "use": "sig",
                "alg": "RS256",
                "kid": kid,
                "n": _b64url(n_bytes),
                "e": _b64url(e_bytes),
            }
        ]
    }
    now = int(time.time())
    hdr = _b64url(json.dumps({"alg": "RS256", "typ": "JWT", "kid": kid}).encode())
    pay = _b64url(
        json.dumps(
            {
                "iss": "https://ubuntubio.org",
                "aud": "openxet",
                "iat": now,
                "nbf": now,
                "exp": now + 3600,
                "jti": "ckpt-1",
                "type": "audit_checkpoint",
                "repo": repo,
                "head_hash": head_hash,
                "entries_checked": entries_checked,
            }
        ).encode()
    )
    signing_input = f"{hdr}.{pay}".encode("ascii")
    sig = subprocess.run(
        ["openssl", "dgst", "-sha256", "-sign", str(pem_file)],
        input=signing_input,
        check=True,
        capture_output=True,
    ).stdout
    return f"{hdr}.{pay}.{_b64url(sig)}", jwks, pem_file


def test_verify_rs256_checkpoint_and_detect_privileged_rewrite(tmp_path: Path) -> None:
    """Detect history rewrite or truncation relative to an externally retained anchor."""
    h1 = hashlib.sha256(b"entry-1").hexdigest()
    h2 = hashlib.sha256(b"entry-2").hexdigest()
    token1, jwks, priv = _sign_test_checkpoint_jwt(
        key_dir=tmp_path,
        repo="waccbip/cholera",
        head_hash=h2,
        entries_checked=2,
    )
    claims1 = verify_rs256_checkpoint_jwt(token1, jwks, expected_repo="waccbip/cholera")
    assert claims1["entries_checked"] == 2
    assert claims1["head_hash"] == h2

    # Tampered JWT signature is rejected
    tampered_tok = token1[:-4] + ("AAAA" if not token1.endswith("AAAA") else "BBBB")
    with pytest.raises(UbuntuBioError, match="Invalid RS256 signature"):
        verify_rs256_checkpoint_jwt(tampered_tok, jwks, expected_repo="waccbip/cholera")

    anchor = [{"repo": "waccbip/cholera", "entries_checked": 2, "head_hash": h2}]

    # 1. Legitimate forward append (entries=3, entry #2 still has hash h2) passes
    h3 = hashlib.sha256(b"entry-3").hexdigest()
    token2, _, _ = _sign_test_checkpoint_jwt(
        key_dir=tmp_path,
        repo="waccbip/cholera",
        head_hash=h3,
        entries_checked=3,
        key_path=priv,
    )
    claims2 = verify_rs256_checkpoint_jwt(token2, jwks, expected_repo="waccbip/cholera")
    clean_log = [
        {"sequence_no": 1, "entry_hash": h1},
        {"sequence_no": 2, "entry_hash": h2},
        {"sequence_no": 3, "entry_hash": h3},
    ]
    assert verify_checkpoint_continuity(claims2, anchor, audit_entries=clean_log) == []

    # 2. Privileged rewrite+append attack (admin alters entry #1, recomputes chain, appends #3,
    #    and re-signs checkpoint with server RSA key) is CAUGHT against external anchor
    forged_h1 = hashlib.sha256(b"tampered-1").hexdigest()
    forged_h2 = hashlib.sha256(b"recomputed-2").hexdigest()
    forged_h3 = hashlib.sha256(b"recomputed-3").hexdigest()
    forged_tok, _, _ = _sign_test_checkpoint_jwt(
        key_dir=tmp_path,
        repo="waccbip/cholera",
        head_hash=forged_h3,
        entries_checked=3,
        key_path=priv,
    )
    forged_claims = verify_rs256_checkpoint_jwt(forged_tok, jwks, expected_repo="waccbip/cholera")
    forged_log = [
        {"sequence_no": 1, "entry_hash": forged_h1},
        {"sequence_no": 2, "entry_hash": forged_h2},
        {"sequence_no": 3, "entry_hash": forged_h3},
    ]
    violations = verify_checkpoint_continuity(forged_claims, anchor, audit_entries=forged_log)
    assert len(violations) == 1
    assert "Historical audit rewrite detected" in violations[0]

    # 3. Chain truncation attack is CAUGHT
    trunc_tok, _, _ = _sign_test_checkpoint_jwt(
        key_dir=tmp_path,
        repo="waccbip/cholera",
        head_hash=h1,
        entries_checked=1,
        key_path=priv,
    )
    trunc_claims = verify_rs256_checkpoint_jwt(trunc_tok, jwks, expected_repo="waccbip/cholera")
    trunc_violations = verify_checkpoint_continuity(trunc_claims, anchor)
    assert any("truncation" in v.lower() for v in trunc_violations)


def test_audit_checkpoint_cli_save_and_verify(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """CLI ubio audit checkpoint saves JSONL external ledger and detects divergence."""
    h1 = hashlib.sha256(b"row-1").hexdigest()
    h2 = hashlib.sha256(b"row-2").hexdigest()
    tok1, jwks, priv = _sign_test_checkpoint_jwt(
        key_dir=tmp_path,
        repo="waccbip/cholera",
        head_hash=h2,
        entries_checked=2,
    )
    audit_rows = [
        {"sequence_no": 1, "entry_hash": h1},
        {"sequence_no": 2, "entry_hash": h2},
    ]

    class _FakeResp:
        def __init__(self, status_code: int, payload: object) -> None:
            self.status_code = status_code
            self._payload = payload
            self.text = json.dumps(payload)

        def json(self) -> object:
            return self._payload

    def _fake_get(
        url: str, headers: dict[str, str] | None = None, timeout: float = 15.0
    ) -> _FakeResp:
        if url.endswith("/.well-known/jwks.json"):
            return _FakeResp(200, jwks)
        if "/audit-log/checkpoint" in url:
            return _FakeResp(
                200,
                {
                    "repo_slug": "waccbip/cholera",
                    "entries_checked": 2,
                    "head_hash": h2,
                    "kid": "test-kid-1",
                    "checkpoint_jwt": tok1,
                },
            )
        if "/audit-log" in url:
            return _FakeResp(200, audit_rows)
        return _FakeResp(404, {})

    monkeypatch.setattr("ubuntu_bio.audit.httpx.get", _fake_get)
    monkeypatch.setenv("UBUNTUBIO_TOKEN", "ubp_test_token")
    monkeypatch.setenv("UBUNTUBIO_HUB_URL", "https://ubuntubio.org")

    ledger = tmp_path / "external_anchors.jsonl"
    runner = CliRunner()
    res_save = runner.invoke(app, ["audit", "checkpoint", "waccbip/cholera", "--save", str(ledger)])
    assert res_save.exit_code == 0, res_save.stdout
    assert "[VERIFIED]" in res_save.stdout
    assert len(load_anchor_ledger(ledger)) == 1

    # Verify against rewritten log fails exit_code=1
    audit_rows[1]["entry_hash"] = hashlib.sha256(b"rewritten-row-2").hexdigest()
    tok2, _, _ = _sign_test_checkpoint_jwt(
        key_dir=tmp_path,
        repo="waccbip/cholera",
        head_hash=hashlib.sha256(b"rewritten-row-2").hexdigest(),
        entries_checked=2,
        key_path=priv,
    )

    def _fake_get_forged(
        url: str, headers: dict[str, str] | None = None, timeout: float = 15.0
    ) -> _FakeResp:
        if "/audit-log/checkpoint" in url:
            return _FakeResp(
                200,
                {
                    "repo_slug": "waccbip/cholera",
                    "entries_checked": 2,
                    "head_hash": hashlib.sha256(b"rewritten-row-2").hexdigest(),
                    "kid": "test-kid-1",
                    "checkpoint_jwt": tok2,
                },
            )
        return _fake_get(url, headers=headers, timeout=timeout)

    monkeypatch.setattr("ubuntu_bio.audit.httpx.get", _fake_get_forged)
    res_fail = runner.invoke(
        app,
        ["audit", "checkpoint", "waccbip/cholera", "--verify-against", str(ledger)],
    )
    assert res_fail.exit_code == 1
    assert "AUDIT FAILED" in res_fail.stdout

    # Direct function call raises on violation
    with pytest.raises(UbuntuBioError, match="Audit rewrite detected"):
        fetch_and_verify_checkpoint(
            "waccbip/cholera",
            verify_against=ledger,
            token="ubp_test_token",
            hub_url="https://ubuntubio.org",
        )
