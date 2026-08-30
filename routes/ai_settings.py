"""Platform AI provider settings routes."""
from flask import Blueprint, jsonify, render_template, request

from modules.ai_config import (
    get_effective_ai_config,
    public_ai_config,
    reset_ai_config,
    save_ai_config,
    test_ai_connection,
)
from routes.auth import require_ai_token


ai_settings_bp = Blueprint("ai_settings", __name__)


@ai_settings_bp.route("/settings/ai")
def ai_settings_page():
    return render_template("ai_settings.html")


@ai_settings_bp.route("/api/settings/ai")
@require_ai_token
def ai_settings_status():
    return jsonify(public_ai_config())


@ai_settings_bp.route("/api/settings/ai", methods=["POST"])
@require_ai_token
def ai_settings_save():
    data = request.get_json(silent=True) or {}
    try:
        values = save_ai_config(
            api_url=data.get("api_url", ""),
            model=data.get("model", ""),
            provider=data.get("provider", "auto"),
            api_key=data.get("api_key"),
            clear_api_key=bool(data.get("clear_api_key", False)),
        )
        return jsonify({"ok": True, "settings": public_ai_config(), "provider": values["provider"]})
    except (ValueError, TypeError) as exc:
        return jsonify({"error": str(exc)}), 400
    except Exception as exc:
        return jsonify({"error": f"保存 AI 配置失败：{exc}"}), 500


@ai_settings_bp.route("/api/settings/ai/reset", methods=["POST"])
@require_ai_token
def ai_settings_reset():
    try:
        reset_ai_config()
        return jsonify({"ok": True, "settings": public_ai_config()})
    except Exception as exc:
        return jsonify({"error": f"恢复默认配置失败：{exc}"}), 500


@ai_settings_bp.route("/api/settings/ai/test", methods=["POST"])
@require_ai_token
def ai_settings_test():
    data = request.get_json(silent=True) or {}
    current = get_effective_ai_config()
    values = {
        "api_url": str(data.get("api_url") or current["api_url"]).strip(),
        "model": str(data.get("model") or current["model"]).strip(),
        "provider": str(data.get("provider") or current["provider"]).strip(),
        # An omitted/blank field means “use the currently saved key”; the key
        # is never included in the response.
        "api_key": str(data.get("api_key") or current["api_key"]).strip(),
    }
    try:
        return jsonify(test_ai_connection(values))
    except (ValueError, TypeError) as exc:
        return jsonify({"ok": False, "error": str(exc)}), 400
    except Exception as exc:
        return jsonify({"ok": False, "error": f"连接测试失败：{exc}"}), 502
