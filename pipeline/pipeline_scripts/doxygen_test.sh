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

#Определяем измененные файлы и записываем их в переменную
if [ "${CI_PIPELINE_SOURCE}" = "merge_request_event" ]; then
  files="$(git diff --name-only origin/${CI_MERGE_REQUEST_SOURCE_BRANCH_NAME} origin/${CI_MERGE_REQUEST_TARGET_BRANCH_NAME}| xargs -I {} basename {} | grep -owP "[^\s]+\.cpp|[^\s]+\.h")"
  echo -e "\033[32m---------------MR---------------"
  echo -e "\033[32morigin/${CI_MERGE_REQUEST_SOURCE_BRANCH_NAME} to origin/${CI_MERGE_REQUEST_TARGET_BRANCH_NAME}"
else
  files="$(git diff-tree --no-commit-id --name-only -r $CI_COMMIT_SHA| xargs -I {} basename {} | grep -owP "[^\s]+\.cpp|[^\s]+\.h")"
  echo -e "\033[32m-------------COMMIT-------------"
fi

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
