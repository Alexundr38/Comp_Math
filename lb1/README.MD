# Comp_Math lb1

--------------------------------------------------------------------------------------------------

-- Install command:

sudo apt-get update

sudo apt-get install build-essential libsdl2-dev libgl1-mesa-dev libglu1-mesa-dev libsdl2-ttf-dev

--------------------------------------------------------------------------------------------------

Run thomas_attractor or use cmake

--------------------------------------------------------------------------------------------------

-- Cmake command:

mkdir build

cd build

cmake ..

make

--------------------------------------------------------------------------------------------------

-- If have troubles with sdl2_ttf.pc (Ubuntu)

1. Create file:

sudo nano /usr/lib/x86_64-linux-gnu/pkgconfig/sdl2_ttf.pc

2. Insert:

prefix=/usr

exec_prefix=${prefix}

libdir=${exec_prefix}/lib/x86_64-linux-gnu

includedir=${prefix}/include

Name: SDL2_ttf

Description: SDL2 TrueType font library

Version: 2.22.0

Requires: sdl2 >= 2.0.8

Libs: -L${libdir} -lSDL2_ttf

Cflags: -I${includedir}/SDL2

3. Save file

4. Reload PKG_CONFIG_PATH:

export PKG_CONFIG_PATH=/usr/lib/x86_64-linux-gnu/pkgconfig:$PKG_CONFIG_PATH

5. If this command output "-lSDL2_ttf" - good result

pkg-config --libs sdl2_ttf

--------------------------------------------------------------------------------------------------

-- Key assigments:

e - use Euler method

m - use MiddlePoint method

r - use RK4 method

t - use Trapezoid method

p - use PredictCorrector 

d - use DOPRI8 method


'>' - increase h

< - decrease h

] - increase b

[ - decrease b


up - rotate up

down - rotate down

left - rotate left

right - rotate right


PgUp - zoom

PgDn - move away



