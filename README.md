# Лабораторная работа: Кубический сплайн

## Сборка проекта

Проект собирается с помощью CMake и MinGW.

1. Убедитесь, что установлены:
   - CMake (версия 3.10 или выше)
   - MinGW-w64 (с компилятором g++)

2. Откройте терминал в корневой папке проекта (там, где лежат `CMakeLists.txt`, `algorithm.cpp`, `accuracy_engine_simple.cpp` и остальные файлы).

3. Выполните команды:
   ```powershell
   mkdir build
   cd build
   cmake .. -G "MinGW Makefiles"
   mingw32-make
   
4. Для запуска и генерации таблиц (в папку билда):
   ```powershell
   ./algorithm_test.exe

5. Для запуска тестов:
   ```powershell
   ./algorithm_test.exe