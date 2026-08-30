"""Runtime AI provider settings.

The process environment remains the default source of configuration.  Values
saved from the platform UI are stored in SQLite and merged at request time so
that changing a provider does not require restarting Flask.
"""
from __future__ import annotations

from typing import Any, Dict

from config import Config


SETTING_KEYS = {
    "api_url": "ai_api_url",
    "api_key": "ai_api_key",
    "model": "ai_model",
    "provider": "ai_provider",
}


def _defaults() -> Dict[str, str]:
    return {
        "api_url": str(Config.AI_API_URL or ""),
        "api_key": str(Config.AI_API_KEY or ""),
        "model": str(Config.AI_MODEL or ""),
        "provider": "auto",
    }


def _read_saved() -> Dict[str, str]:
    """Read saved settings, returning an empty mapping when DB is unavailable."""
    try:
        from database import get_conn

        conn = get_conn()
        try:
            rows = conn.execute(
                "SELECT key, value FROM platform_settings WHERE key IN (?, ?, ?, ?)",
                tuple(SETTING_KEYS.values()),
            ).fetchall()
        finally:
            conn.close()
        reverse = {value: key for key, value in SETTING_KEYS.items()}
        return {reverse[row["key"]]: (row["value"] or "") for row in rows}
    except Exception:
        # Config must still work during early startup and in isolated scripts
        # that do not initialize the application database.
        return {}


def get_effective_ai_config() -> Dict[str, str]:
    """Return the currently effective, non-secret AI configuration."""
    values = _defaults()
    values.update(_read_saved())
    return values


def save_ai_config(*, api_url: str, model: str, provider: str = "auto",
                   api_key: str | None = None, clear_api_key: bool = False) -> Dict[str, str]:
    """Persist UI settings and return the effective values.

    A blank key keeps the existing key by default, which lets the UI show a
    masked key without requiring the user to paste it again on every edit.
    ``clear_api_key`` explicitly removes it.
    """
    api_url = str(api_url or "").strip()
    model = str(model or "").strip()
    provider = str(provider or "auto").strip().lower()
    if not api_url:
        raise ValueError("API 地址不能为空")
    if not model:
        raise ValueError("模型名不能为空")
    if provider not in {"auto", "openai", "anthropic"}:
        raise ValueError("协议必须是 auto、openai 或 anthropic")
    if not (api_url.startswith("http://") or api_url.startswith("https://")):
        raise ValueError("API 地址必须以 http:// 或 https:// 开头")

    current = get_effective_ai_config()
    if clear_api_key:
        effective_key = ""
    elif api_key is None or not str(api_key).strip():
        effective_key = current["api_key"]
    else:
        effective_key = str(api_key).strip()

    values = {
        "api_url": api_url,
        "api_key": effective_key,
        "model": model,
        "provider": provider,
    }
    from database import get_conn

    conn = get_conn()
    try:
        for name, value in values.items():
            conn.execute(
                "INSERT INTO platform_settings (key, value, updated_at) VALUES (?, ?, CURRENT_TIMESTAMP) "
                "ON CONFLICT(key) DO UPDATE SET value=excluded.value, updated_at=CURRENT_TIMESTAMP",
                (SETTING_KEYS[name], value),
            )
        conn.commit()
    finally:
        conn.close()
    return values


def reset_ai_config() -> Dict[str, str]:
    """Remove UI overrides and fall back to environment/config defaults."""
    from database import get_conn

    conn = get_conn()
    try:
        conn.execute(
            "DELETE FROM platform_settings WHERE key IN (?, ?, ?, ?)",
            tuple(SETTING_KEYS.values()),
        )
        conn.commit()
    finally:
        conn.close()
    return get_effective_ai_config()


def mask_api_key(value: str | None) -> str:
    value = str(value or "")
    if not value:
        return "未配置"
    if len(value) <= 6:
        return "••••••"
    return f"{value[:3]}{'•' * max(4, len(value) - 6)}{value[-3:]}"


def public_ai_config() -> Dict[str, Any]:
    values = get_effective_ai_config()
    return {
        "api_url": values["api_url"],
        "model": values["model"],
        "provider": values["provider"],
        "configured": bool(values["api_key"]),
        "api_key_masked": mask_api_key(values["api_key"]),
        "has_saved_override": bool(_read_saved()),
        "requires_local_token": bool(Config.AI_API_TOKEN),
    }


def test_ai_connection(values: Dict[str, str]) -> Dict[str, Any]:
    """Perform a small provider health request without changing saved config."""
    import requests

    url = str(values.get("api_url") or "").strip()
    key = str(values.get("api_key") or "").strip()
    model = str(values.get("model") or "").strip()
    provider = str(values.get("provider") or "auto").strip().lower()
    if not key:
        raise ValueError("请先填写 API Key（测试连接不会保存配置）")
    if not url or not model:
        raise ValueError("API 地址和模型名不能为空")
    is_anthropic = provider == "anthropic" or (
        provider == "auto" and ("anthropic" in url.lower() or "claude" in url.lower())
    )
    if is_anthropic:
        endpoint = url.rstrip("/")
        if not endpoint.endswith(("/messages", "/v1/messages")):
            endpoint += "/messages" if endpoint.endswith("/v1") else "/v1/messages"
        response = requests.post(
            endpoint,
            headers={
                "Content-Type": "application/json",
                "x-api-key": key,
                "Authorization": f"Bearer {key}",
                "anthropic-version": "2023-06-01",
            },
            json={"model": model, "max_tokens": 1, "messages": [{"role": "user", "content": "ping"}]},
            timeout=10.0,
        )
    else:
        endpoint = url.rstrip("/")
        for suffix in ("/chat/completions", "/v1/chat/completions"):
            if endpoint.endswith(suffix):
                endpoint = endpoint[:-len(suffix)]
                break
        response = requests.get(
            endpoint + ("/models" if endpoint.endswith("/v1") else "/v1/models"),
            headers={"Authorization": f"Bearer {key}", "Content-Type": "application/json"},
            timeout=10.0,
        )
    if response.status_code >= 400:
        detail = response.text[:240]
        raise RuntimeError(f"服务返回 HTTP {response.status_code}: {detail}")
    return {"ok": True, "provider": "anthropic" if is_anthropic else "openai", "status_code": response.status_code}
