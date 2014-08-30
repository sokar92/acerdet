CC= g++
CFLAGS= -g3 -O2 -fPIC -Wall

HEPMC_CONVERTER= HepMC_InputConverter
PYTHIA_PROVIDER= Pythia8_ParticleDataProviderFactory
ROOT_HISTOGRAM= Root_HistogramManager

OBJECTS= $(HEPMC_CONVERTER).o \
	$(PYTHIA_PROVIDER).o \
	$(ROOT_HISTOGRAM).o

ROOTINCDIR := $(shell root-config --incdir)

all: prepare $(OBJECTS) export_headers

$(HEPMC_CONVERTER).o : $(HEPMC_CONVERTER).cpp
	$(CC) $(CFLAGS) -c $< -o $@ -I$(HEPMCLOCATION)/include

$(PYTHIA_PROVIDER).o : $(PYTHIA_PROVIDER).cpp
	$(CC) $(CFLAGS) -c $< -o $@ -I$(PYTHIALOCATION)/include

$(ROOT_HISTOGRAM).o : $(ROOT_HISTOGRAM).cpp
	$(CC) $(CFLAGS) -c $< -o $@ -I$(ROOTINCDIR)

INPUT_HEADERS_PATH= ../include/external

prepare:
	mkdir $(INPUT_HEADERS_PATH)

export_headers:
	cp -f *.h -t $(INPUT_HEADERS_PATH)

clean:
	find . -regex '.*\.\(~\|o\|so\)' | xargs rm -rf
