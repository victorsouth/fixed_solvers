#!/bin/bash

# define paths
PWD_STORE=$PWD

# Проверяем, передан ли аргумент (GCC, MinGW или Clang)
if [ -z "$1" ]; then
    echo "Usage: $0 [GCC|MinGW|Clang]"
    exit 1
fi

COMPILER=$1

case $COMPILER in
    GCC)
        echo "Building with GCC"
        LOG_DIR="$PWD/pipeline/pipeline_result/build_gcc_result"
		LIB_DIR="$PWD/Libs/GCC"
		INSTALL_DIR="$PWD/Libs/GCC"
        ;;
    MinGW)
        echo "Building with MinGW in MSYS2..."
        LOG_DIR="$PWD/pipeline/pipeline_result/build_mingw_result"
		LIB_DIR="$PWD/Libs/MinGW"
		INSTALL_DIR="$PWD/Libs/MinGW"
        ;;
    Clang)
        echo "Building with Clang in MSYS2..."
        LOG_DIR="$PWD/pipeline/pipeline_result/build_clang_result"
		LIB_DIR="$PWD/Libs/Clang"
		INSTALL_DIR="$PWD/Libs/Clang"
        ;;
    *)
        echo "Invalid compiler! Use GCC, MinGW, or Clang."
        exit 1
        ;;
esac

# ----------------------------------------------------------------------------------------------------------------------------------------------
# Сборка и установка
# ----------------------------------------------------------------------------------------------------------------------------------------------
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="-w -Wno-#pragma-messages" -DCMAKE_INSTALL_PREFIX="$INSTALL_DIR" > ${LOG_DIR}/fixed_solvers_build.log 2>&1
cmake --build . -j$(nproc) >> ${LOG_DIR}/fixed_solvers_build.log 2>&1
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
cmake --install . >> ${LOG_DIR}/fixed_solvers_build.log 2>&1
if [ $? -ne 0 ]; then
    echo "--------------- Result ---------------"
	echo "Ошибка: установка fixed_solvers [$FIXED_SOLVERS_BRANCH] не удалась"
	echo "--------------------------------------"
    exit 1
fi
echo "fixed_solvers [$FIXED_SOLVERS_BRANCH] успешно установлен"
cd $PWD_STORE