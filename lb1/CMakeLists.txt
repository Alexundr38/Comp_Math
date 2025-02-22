cmake_minimum_required(VERSION 3.10)

project(Thomas_Attractor)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR})

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(SDL2 REQUIRED)
find_package(OpenGL REQUIRED)
find_package(PkgConfig REQUIRED)
pkg_check_modules(SDL2_TTF REQUIRED sdl2_ttf)

add_executable(thomas_attractor main.cpp)

target_link_libraries(thomas_attractor
    ${SDL2_LIBRARIES}
    ${OPENGL_LIBRARIES}
    ${SDL2_TTF_LIBRARIES}
)

target_include_directories(thomas_attractor PRIVATE
    ${SDL2_INCLUDE_DIRS}
    ${OPENGL_INCLUDE_DIRS}
    ${SDL2_TTF_INCLUDE_DIRS}
)
