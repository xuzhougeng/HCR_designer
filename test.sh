#!/bin/bash
# Run the whole unit-test suite from the project root.
# Uses HCR_PYTHON when set, otherwise the project virtualenv, otherwise python3.

# 设置颜色代码
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# 打印带颜色的消息
print_message() {
    color=$1
    message=$2
    echo -e "${color}${message}${NC}"
}

cd "$(dirname "$0")" || exit 1

if [ -n "$HCR_PYTHON" ]; then
    PYTHON="$HCR_PYTHON"
elif [ -x ".venv/bin/python" ]; then
    PYTHON=".venv/bin/python"
else
    PYTHON="python3"
fi

# 检查Python是否安装
if ! command -v "$PYTHON" &> /dev/null; then
    print_message "$RED" "错误: 未找到 Python ($PYTHON)"
    exit 1
fi

# 检查必要的Python包是否安装
print_message "$YELLOW" "检查依赖项..."
if ! "$PYTHON" -c "import Bio, flask" 2>/dev/null; then
    print_message "$RED" "错误: 未找到 Biopython 或 Flask"
    print_message "$YELLOW" "请运行: $PYTHON -m pip install -r requirements.txt"
    exit 1
fi

# 运行测试
print_message "$YELLOW" "开始运行测试..."
echo "----------------------------------------"

"$PYTHON" -m unittest discover -s tests "$@"
test_result=$?

if command -v node &> /dev/null; then
    node --test tests/test_result_evidence.js || test_result=1
else
    print_message "$YELLOW" "未找到 node，跳过 tests/test_result_evidence.js"
fi

echo "----------------------------------------"

# 根据测试结果显示相应消息
if [ $test_result -eq 0 ]; then
    print_message "$GREEN" "✓ 所有测试通过！"
else
    print_message "$RED" "✗ 测试失败！"
fi

exit $test_result
