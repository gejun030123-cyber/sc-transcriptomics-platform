#!/usr/bin/env bash
# Manage the web platform without interrupting the independent NAT123 tunnel.

set -Eeuo pipefail

SERVICE_NAME="${PLATFORM_SERVICE_NAME:-sc-transcriptomics-platform.service}"
HEALTH_URL="${PLATFORM_HEALTH_URL:-http://127.0.0.1:${PORT:-5000}/login}"
ACTION="${1:-restart}"

run_privileged() {
    if (( EUID == 0 )); then
        "$@"
    elif command -v sudo >/dev/null 2>&1; then
        sudo "$@"
    else
        echo "需要 root 权限执行：$*" >&2
        exit 1
    fi
}

show_status() {
    systemctl status "$SERVICE_NAME" --no-pager -l || true
}

http_is_ready() {
    command -v curl >/dev/null 2>&1 && \
        curl --fail --silent --show-error --output /dev/null \
            --max-time 2 "$HEALTH_URL" 2>/dev/null
}

refuse_unmanaged_listener() {
    if systemctl is-active --quiet "$SERVICE_NAME"; then
        return 0
    fi
    if http_is_ready; then
        echo "检测到 $HEALTH_URL 已被非 systemd 平台进程占用。" >&2
        echo "为避免误杀正在运行的分析，本脚本不会自动终止它。" >&2
        echo "确认没有分析任务后，先停止旧的 python app.py 进程，再重新执行本命令。" >&2
        return 1
    fi
}

wait_until_ready() {
    local attempt
    if ! command -v curl >/dev/null 2>&1; then
        echo "平台服务已启动；未找到 curl，跳过 HTTP 就绪检查。"
        return 0
    fi

    for (( attempt = 1; attempt <= 60; attempt++ )); do
        if systemctl is-active --quiet "$SERVICE_NAME" && http_is_ready; then
            echo "平台已就绪：$HEALTH_URL"
            echo "nat123 服务未重启，原外网地址继续使用。"
            return 0
        fi
        if systemctl is-failed --quiet "$SERVICE_NAME"; then
            echo "平台启动失败，当前服务状态：" >&2
            show_status >&2
            return 1
        fi
        sleep 1
    done

    echo "等待平台就绪超时（$HEALTH_URL），当前服务状态：" >&2
    show_status >&2
    return 1
}

case "$ACTION" in
    start)
        refuse_unmanaged_listener
        run_privileged systemctl start "$SERVICE_NAME"
        wait_until_ready
        ;;
    restart|reload)
        refuse_unmanaged_listener
        run_privileged systemctl restart "$SERVICE_NAME"
        wait_until_ready
        ;;
    stop)
        run_privileged systemctl stop "$SERVICE_NAME"
        echo "平台已停止；nat123 服务未停止。"
        ;;
    status)
        show_status
        ;;
    logs)
        run_privileged journalctl -u "$SERVICE_NAME" -n 100 -f
        ;;
    help|-h|--help)
        echo "用法：$0 {restart|start|stop|status|logs}"
        echo "默认动作：restart；所有动作只管理平台，不管理 nat123。"
        ;;
    *)
        echo "未知动作：$ACTION" >&2
        echo "用法：$0 {restart|start|stop|status|logs}" >&2
        exit 2
        ;;
esac
