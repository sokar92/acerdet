CC= g++
CFLAGS= -g3 -O2 -fPIC -Wall

OBJECTS= Histogram.o Particle.o Functions.o 

all: prepare $(OBJECTS) export_include

%.o : %.cpp
	$(CC) $(CFLAGS) -c $< -o $@ -lStd

CORE_HEADERS_PATH= ../include/core
prepare:
	mkdir $(CORE_HEADERS_PATH)

export_include: 
	cp -f *.h -t $(CORE_HEADERS_PATH)

clean:
	find . -regex '.*\.\(~\|o\|so\)' | xargs rm -rf
