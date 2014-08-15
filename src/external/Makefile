CC= g++
CFLAGS= -g3 -O2 -fPIC -Wall

OBJECTS= HepMC_InputConverter.o

all: prepare $(OBJECTS) export_headers

%.o : %.cpp
	$(CC) $(CFLAGS) -c $< -o $@ -I../$(HEPMCLOCATION)include

INPUT_HEADERS_PATH= ../include/io

prepare:
	mkdir $(INPUT_HEADERS_PATH)

export_headers:
	cp -f *.h -t $(INPUT_HEADERS_PATH)

clean:
	find . -regex '.*\.\(~\|o\|so\)' | xargs rm -rf
