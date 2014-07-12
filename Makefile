CC= g++
CFLAGS= -g3 -O2 -fPIC -Wall
DEMO= demo.exe demoCreateConfig.exe
HEPMCLIBS= -L$(HEPMCLOCATION)lib -lHepMC
PYTHIALIBS= -L$(PYTHIALOCATION)lib -lpythia8 -llhapdfdummy -lpythia8tohepmc

all: $(DEMO)

%.o : %.cpp
	$(CC) $(CFLAGS) -c $< -o $@ -I$(PYTHIALOCATION)include -I$(HEPMCLOCATION)include

%.exe : %.o
	$(CC) $< -o $@ -lAcerDET $(PYTHIALIBS) $(HEPMCLIBS) -Wl,-rpath,$(PYTHIALOCATION)lib -Wl,-rpath,$(HEPMCLOCATION)lib

clean:
	rm -f *~ *.o *.exe
