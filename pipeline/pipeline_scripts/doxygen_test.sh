#!/bin/bash
# Скрипт предназначен для поиска измененных файлов, идущих на MR, и прогонки их через doxygen.

# Забираем изменения из репозитория
git fetch

#Получаем абсолютный путь к директории, в которой находится этот скрипт
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Забираем файл конфигурации и переносим в корневую папку
cp "${SCRIPT_DIR}/Doxyfile" "${SCRIPT_DIR}/../../Doxyfile"

#Получаем путь к файлу конфигурации doxygen
Doxyfile="${SCRIPT_DIR}/../../Doxyfile"

for file in $(find . -type f \(-name "*.cpp" -o -name "*.h" \)); do
	echo $file >> files
done

#Выходим из скрипта, если файл с изменениями пуст
if [[ -z "$files" ]]; then
	echo -e "\033[32m--------No-files-to-check--------"
	exit 0
fi

echo "--------------FILES--------------"
echo -e "\033[32m$files"

# Записываем в Doxyfile перечень измененных файлов в раздел FILE_PATTERNS
first_line=true
for file in $files
do
  if $first_line; then
    echo "FILE_PATTERNS = $file \\" >> ${Doxyfile}
    first_line=false
  else
    echo " $file \\" >> ${Doxyfile}
  fi
done

# Вызываем команду doxygen
doxygen ${Doxyfile} > doxygen_out.log 2>&1
cp "${SCRIPT_DIR}/../../doxygen_out.log" "${SCRIPT_DIR}/../pipeline_result/doxygen_out"
cp "${SCRIPT_DIR}/../../doxygen_test_out" "${SCRIPT_DIR}/../pipeline_result/doxygen_test_out"

# Записываем результат в doxygen_test_out.txt и выводим его, если не пустой

if [ -s ${SCRIPT_DIR}/../../doxygen_test_out ]; then
	echo -e "\033[33m----------Warnings----------"
	tac ${SCRIPT_DIR}/../../doxygen_test_out | sort
    exit 1;
else
	echo -e "\033[32m--------No-comment-warnings---------"
	exit 0
fi
