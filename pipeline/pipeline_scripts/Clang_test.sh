#!/bin/bash
# Скрипт для статического анализа всей кодовой базы
git fetch
#Получаем абсолютный путь к директории, в которой находится этот скрипт
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

for file in $(find . -type f \(-name "*.cpp" -o -name "*.h" \)); do
	echo $file
	clang-tidy $file -checks='modernize-use-default-member-init,modernize-use-default,cppcoreguidelines-virtual-class-destructor,modernize-use-using,readability-convert-member-functions-to-static' >> clang_full_log.txt 2>&1
done

echo "============full:=============="
cat "clang_full_log.txt" | grep "modernize-use-default-member-init\|modernize-use-default\|cppcoreguidelines-virtual-class-destructor\|modernize-use-using\|readability-convert-member-functions-to-static"| egrep -v "clang-diagnostic" >> clang_log.txt

cp "clang_log.txt" "${SCRIPT_DIR}/../pipeline_result/clang_log.txt"

if [ -s clang_log.txt ]; then
	echo -e "\033[33m----------Warnings----------"
	tac clang_log.txt
    exit 1;
else
	echo -e "\033[32m--------No-warnings---------"
	exit 0
fi