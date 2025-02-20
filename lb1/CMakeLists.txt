# Минимальная требуемая версия CMake
cmake_minimum_required(VERSION 3.10)

# Название проекта
project(OpenGL_SDL_Example)

# Установка стандарта C++ (например, C++14)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Поиск необходимых библиотек
find_package(SDL2 REQUIRED)
find_package(OpenGL REQUIRED)
find_package(PkgConfig REQUIRED)
pkg_check_modules(SDL2_TTF REQUIRED sdl2_ttf)

# Добавление исполняемого файла
add_executable(opengl_sdl_example main.cpp)

# Линковка библиотек
target_link_libraries(opengl_sdl_example
    ${SDL2_LIBRARIES}
    ${OPENGL_LIBRARIES}
    ${SDL2_TTF_LIBRARIES}
)

# Включение директорий с заголовочными файлами
target_include_directories(opengl_sdl_example PRIVATE
    ${SDL2_INCLUDE_DIRS}
    ${OPENGL_INCLUDE_DIRS}
    ${SDL2_TTF_INCLUDE_DIRS}
)
