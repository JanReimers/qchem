//
//  In order to set up the std headers for modules you need run a command in build directory like this:
//      cd build/Release
//      g++-15 -v -c -fmodules-ts -x c++-system-header  -std=c++20 iostream string vector
//  Add more header as required.
//
// To help VSCode understand C++20 module imports, ensure your project uses a recent compiler (like GCC 13+ or Clang 16+), 
// and configure your build system (CMake, tasks.json, etc.) to use -std=c++20 and -fmodules-ts. 
// Also, install the latest C/C++ extension for VSCode. 
// For IntelliSense, set "C_Cpp.default.cppStandard": "c++20" in your settings.json.
// Note: Full module support in VSCode IntelliSense is still experimental and may require future updates.

import <iostream>;

int main() {
    std::cout << "Hello, World! This is a C++20 program." << std::endl;
    return 0;
}