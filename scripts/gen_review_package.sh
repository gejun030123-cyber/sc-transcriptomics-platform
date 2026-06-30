#!/bin/bash
cd "$(dirname "$0")/.."
OUTPUT="review_package.md"

echo "# 生信分析平台 — 代码审查包" > "$OUTPUT"
echo "" >> "$OUTPUT"
echo "生成时间: $(date '+%Y-%m-%d %H:%M')" >> "$OUTPUT"
echo "" >> "$OUTPUT"

echo "## 项目结构" >> "$OUTPUT"
echo '```' >> "$OUTPUT"
find . -name "*.py" -not -path "./.git/*" -not -path "./__pycache__/*" -not -path "./tests/*" | sort >> "$OUTPUT"
echo '```' >> "$OUTPUT"
echo "" >> "$OUTPUT"

echo "## 最近 10 个 commit" >> "$OUTPUT"
echo '```' >> "$OUTPUT"
git log --oneline -10 >> "$OUTPUT"
echo '```' >> "$OUTPUT"
echo "" >> "$OUTPUT"

if [ $# -gt 0 ]; then
    for f in "$@"; do
        if [ -f "$f" ]; then
            echo "## 文件: $f" >> "$OUTPUT"
            echo '```python' >> "$OUTPUT"
            cat "$f" >> "$OUTPUT"
            echo '```' >> "$OUTPUT"
            echo "" >> "$OUTPUT"
        fi
    done
else
    echo "## 用法" >> "$OUTPUT"
    echo 'bash scripts/gen_review_package.sh modules/bulk_deg.py modules/ai_adapter.py' >> "$OUTPUT"
fi

echo "✅ 已生成: $OUTPUT ($(wc -l < "$OUTPUT") 行)"
