#!/bin/bash

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

# 检查Python是否安装
if ! command -v python3 &> /dev/null; then
    print_message "$RED" "错误: 未找到 Python3"
    exit 1
fi

# 检查必要的Python包是否安装
print_message "$YELLOW" "检查依赖项..."
python3 -c "import Bio" 2>/dev/null
if [ $? -ne 0 ]; then
    print_message "$RED" "错误: 未找到 Biopython 包"
    print_message "$YELLOW" "请运行: pip install biopython"
    exit 1
fi

# 检查测试文件是否存在
if [ ! -f "tests/test_tcr2.py" ]; then
    print_message "$RED" "错误: 未找到测试文件 tests/test_tcr2.py"
    exit 1
fi

# 检查源代码文件是否存在
if [ ! -f "scripts/tcr2.py" ]; then
    print_message "$RED" "错误: 未找到源代码文件 scripts/tcr2.py"
    exit 1
fi

# 运行测试
print_message "$YELLOW" "开始运行测试..."
echo "----------------------------------------"

# 使用 -v 参数运行测试以获得详细输出
python3 -m unittest tests/test_tcr2.py -v

# 获取测试结果
test_result=$?

echo "----------------------------------------"

# 根据测试结果显示相应消息
if [ $test_result -eq 0 ]; then
    print_message "$GREEN" "✓ 所有测试通过！"
else
    print_message "$RED" "✗ 测试失败！"
fi

exit $test_result
