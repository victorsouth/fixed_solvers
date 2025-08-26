Подготовка Debug

Должны быть предсобраны и заинстоллены с помощью CMAKE:
- Eigen Debug
- Gtest Debug

При всех собранных зависимостях используем типовую команду (запускать из fixed_solvers/msvc_cmake)

cmake .. -G "Visual Studio 17 2022" -DCMAKE_PREFIX_PATH="C:/Program Files (x86)/googletest-distribution"  -DCMAKE_MSVC_RUNTIME_LIBRARY=MultiThreadedDebugDLL -DCMAKE_BUILD_TYPE=Debug

В папке fixed_solvers/msvc_cmake будет создано решение MSVC (.sln). Его открываем, работаем с ним как обычно.

Для установки с целью использования в другом зависимом проекте (pde_solvers и др.):

cmake --build . --config Debug -j 8
cmake --install . --config Debug



Подготовка Release

Должны быть предсобраны и заинстоллены с помощью CMAKE:
- Eigen Release
- Gtest Release

При всех собранных зависимостях используем типовую команду (запускать из fixed_solvers/msvc_cmake)

cmake .. -G "Visual Studio 17 2022" -DCMAKE_PREFIX_PATH="C:/Program Files (x86)/googletest-distribution"  -DCMAKE_MSVC_RUNTIME_LIBRARY=MultiThreadedDLL -DCMAKE_BUILD_TYPE=Release

В папке fixed_solvers/msvc_cmake будет создано решение MSVC (.sln). Его открываем, работаем с ним как обычно.

Для установки с целью использования в другом зависимом проекте (pde_solvers и др.):

cmake --build . --config Release -j 8
cmake --install . --config Release


Сборка зависимостей (команды выполнять в предварительно созданной папке в директории библиотеки с файлом CMakeLists.txt)

-Gtest

DEBUG

cmake  .. -G "Visual Studio 17 2022" -DINSTALL_GTEST=ON -DCMAKE_DEBUG_POSTFIX="d" -DCMAKE_BUILD_TYPE=Debug -Dgtest_force_shared_crt=on
cmake --build . --config Debug -j 8
cmake --install . --config Debug
	

RELEASE
	
cmake  .. -G "Visual Studio 17 2022" -DINSTALL_GTEST=ON -DCMAKE_BUILD_TYPE=Release -Dgtest_force_shared_crt=on
cmake --build . --config Release -j 8
cmake --install . --config Release

-Eigen

cmake  .. -G "Visual Studio 17 2022"
cmake --build . 
cmake --install .



