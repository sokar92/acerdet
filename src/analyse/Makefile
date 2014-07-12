CC= g++
CFLAGS= -g3 -O2 -fPIC -Wall
HEPMCLIBS= -L$(HEPMCLOCATION)lib -lHepMC

OBJECTS= BJet.o \
	Calibration.o \
	Cell.o \
	CJet.o \
	Cluster.o \
	Electron.o \
	Jet.o \
	Mis.o \
	Muon.o \
	Photon.o \
	Tau.o

all: prepare $(OBJECTS) export_analyse

%.o : %.cpp
	$(CC) $(CFLAGS) -c $< -o $@ -lStd

ANALYSE_HEADERS_PATH= ../include/analyse
prepare:
	mkdir $(ANALYSE_HEADERS_PATH)

export_analyse: 
	cp -f *.h -t $(ANALYSE_HEADERS_PATH)

clean:
	find . -regex '.*\.\(~\|o\|so\)' | xargs rm -rf
