#!/bin/sh

#Получаем абсолютный путь к директории, в которой находится этот скрипт
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

#Создаем директорию для файлов журнала
LOG_DIR="${SCRIPT_DIR}/../pipeline_result/build_gcc_result"
mkdir -p ${LOG_DIR}

#Сборка и установка fixed_solvers
mkdir -p build
cd build 
#Собираем fixed_solvers
cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="-w -Wno-#pragma-messages" -DCMAKE_INSTALL_PREFIX=../Libs/GCC > ${LOG_DIR}/fixed_solvers_build.log 2>&1
make -j$(nproc) >> ${LOG_DIR}/fixed_solvers_build.log 2>&1
MAKE_RESULT=$?
# Выводим ошибки из файла fixed_solvers_build.log, содержащие 'error:' или 'fatal error:'
grep -E "error:|fatal error:|Error" ${LOG_DIR}/fixed_solvers_build.log
if [ $MAKE_RESULT -ne 0 ]; then
    echo "--------------- Result ---------------"
	echo "Ошибка: сборка fixed_solvers [$FIXED_SOLVERS_BRANCH] не удалась"
    echo "--------------------------------------"
	exit 1
fi
echo "fixed_solvers [$FIXED_SOLVERS_BRANCH] успешно собран"
#Устанавливаем fixed_solvers
make install >> ${LOG_DIR}/fixed_solvers_build.log 2>&1
if [ $? -ne 0 ]; then
    echo "--------------- Result ---------------"
	echo "Ошибка: установка fixed_solvers [$FIXED_SOLVERS_BRANCH] не удалась"
	echo "--------------------------------------"
    exit 1
fi
echo "fixed_solvers [$FIXED_SOLVERS_BRANCH] успешно установлен"