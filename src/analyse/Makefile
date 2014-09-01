RELATIVE_PATH= ../..
include $(RELATIVE_PATH)/common.inc

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

all: $(OBJECTS)

%.o : %.cpp
	$(CC) $(CFLAGS) -c $< -o $@ -lStd

clean:
	find . -regex '.*\.\(~\|o\|so\)' | xargs rm -rf
