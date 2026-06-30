#!/bin/bash
# =============================================================
# 代码改进 Loop 执行器
# =============================================================
# 用法:
#   bash scripts/run-loop.sh              # 默认最多 20 次迭代
#   bash scripts/run-loop.sh 5            # 最多 5 次迭代
#   bash scripts/run-loop.sh 10 --dry-run # 只打印，不执行
#
# 每次迭代：扫描一个模块 → 生成测试 → 运行验证 → 更新 STATE.md
# =============================================================

set -euo pipefail

MAX_ITER=${1:-20}
DRY_RUN="${2:-}"
STATE_FILE="STATE.md"
PROMPT_FILE="PROMPT.md"
LOG_DIR="logs/loop"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
LOG_FILE="${LOG_DIR}/loop_${TIMESTAMP}.log"

mkdir -p "$LOG_DIR"

log() {
    local msg="[$(date '+%Y-%m-%d %H:%M:%S')] $*"
    echo "$msg"
    echo "$msg" >> "$LOG_FILE"
}

count_remaining() {
    # 统计 STATE.md 中待测模块数量
    grep -c '⬜ 待测' "$STATE_FILE" 2>/dev/null || echo "0"
}

# ---- 主循环 ----
log "========================================="
log "代码改进 Loop 启动"
log "最大迭代: $MAX_ITER"
log "待测模块: $(count_remaining)"
log "日志文件: $LOG_FILE"
log "========================================="

for i in $(seq 1 "$MAX_ITER"); do
    log ""
    log "=== 迭代 $i / $MAX_ITER ==="

    remaining=$(count_remaining)
    if [ "$remaining" -eq 0 ]; then
        log "🎉 所有模块已测试完成！"
        break
    fi
    log "剩余待测模块: $remaining"

    if [ "$DRY_RUN" = "--dry-run" ]; then
        log "[DRY-RUN] 跳过实际执行"
        # 显示下一个待测模块
        grep '⬜ 待测' "$STATE_FILE" | head -1
    else
        log "执行测试生成循环..."
        if claude-code --print "$(cat "$PROMPT_FILE")

当前迭代: $i / $MAX_ITER，剩余 $remaining 个模块。
请执行 Discover 阶段，找到 STATE.md 中下一个优先级最高的待测模块，开始测试生成。" 2>&1 | tee -a "$LOG_FILE"; then
            log "迭代 $i 完成"
        else
            log "⚠️ 迭代 $i 出错（退出码 $?），继续..."
        fi
    fi

    log "等待 3 秒..."
    sleep 3
done

log ""
log "========================================="
log "Loop 结束，共执行 $i 次迭代"
log "日志: $LOG_FILE"
log "========================================="
