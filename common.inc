CC= g++
CFLAGS= -g3 -std=c++0x -O2 -fPIC -Wall -Wno-sign-compare
FASTJETINCLUDE= -I/usr/local/include
FASTJETLIBS= -Wl,-rpath,/usr/local/lib -lm -L/usr/local/lib -lfastjettools -lfastjet
