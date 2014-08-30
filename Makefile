CC= g++
CFLAGS= -g3 -O2 -fPIC -Wall

HEPMCLIBS   = -L$(HEPMCLOCATION)/lib -lHepMC -Wl,-rpath,$(HEPMCLOCATION)/lib
PYTHIALIBS  = -L$(PYTHIALOCATION)/lib -lpythia8 -llhapdfdummy -lpythia8tohepmc -Wl,-rpath,$(PYTHIALOCATION)/lib
ROOTLIBS   := $(shell root-config --libs)

LIBS     = src/libAcerDET.so libHistoManager/libLibHistoManager.so
LIBSLINK = -Lsrc -lAcerDET -Wl,-rpath,src \
           -LlibHistoManager -lLibHistoManager -Wl,-rpath,libHistoManager

all: demo.exe demoCreateConfig.exe

demo.exe: $(LIBS) demo.o
	$(CC) $(CFLAGS) demo.o -o $@ $(LIBSLINK) $(PYTHIALIBS) $(HEPMCLIBS) $(ROOTLIBS)

demoCreateConfig.exe: $(LIBS) demoCreateConfig.o
	$(CC) $(CFLAGS) demoCreateConfig.o -o $@ $(LIBSLINK) $(PYTHIALIBS) $(HEPMCLIBS) $(ROOTLIBS)

src/libAcerDET.so:
	make -C src

libHistoManager/libLibHistoManager.so:
	make -C libHistoManager

%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@ -I$(PYTHIALOCATION)/include -I$(HEPMCLOCATION)/include -Isrc

clean:
	rm -f *~ *.o *.exe *.a

Clean: clean
	make clean -C src
	make clean -C libHistoManager
