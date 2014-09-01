RELATIVE_PATH= ../..
include $(RELATIVE_PATH)/common.inc

OBJECTS= Functions.o \
	CellData.o \
	ClusterData.o \
	JetData.o \
	PartData.o \
	Particle.o \
	ParticleDataProvider.o

all: $(OBJECTS)

%.o : %.cpp
	$(CC) $(CFLAGS) -c $< -o $@ -lStd

clean:
	find . -regex '.*\.\(~\|o\|so\)' | xargs rm -rf
