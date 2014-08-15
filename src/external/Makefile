CC= g++
CFLAGS= -g3 -O2 -fPIC -Wall

OBJECTS= HepMC_InputConverter.o \
	Pythia8_ParticleDataProviderFactory.o \
	Root_HistogramManager.o

all: prepare $(OBJECTS) export_headers

%.o : %.cpp
	$(CC) $(CFLAGS) -c $< -o $@ -I../$(HEPMCLOCATION)include -I$(ROOTSYS)/include

INPUT_HEADERS_PATH= ../include/external

prepare:
	mkdir $(INPUT_HEADERS_PATH)

export_headers:
	cp -f *.h -t $(INPUT_HEADERS_PATH)

clean:
	find . -regex '.*\.\(~\|o\|so\)' | xargs rm -rf
