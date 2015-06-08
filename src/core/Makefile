RELATIVE_PATH= ../..
include $(RELATIVE_PATH)/common.inc

OBJECTS= Smearing.o \
	CellData.o \
	ClusterData.o \
	JetData.o \
	PartData.o \
	MisData.o \
	Particle.o \
	ParticleDataProvider.o \
	Functions.o \
	Kinematics.o

all: $(OBJECTS)

%.o : %.cpp
	$(CC) $(CFLAGS) -c $< -o $@ -lStd

clean:
	find . -regex '.*\.\(~\|o\|so\)' | xargs rm -rf
