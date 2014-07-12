CC= g++
CFLAGS= -g3 -O2 -fPIC -Wall

OBJECTS= Configuration.o 

all: prepare $(OBJECTS) export_include

%.o : %.cpp
	$(CC) $(CFLAGS) -c $< -o $@ -lStd

CONF_HEADERS_PATH= ../include/conf
prepare:
	mkdir $(CONF_HEADERS_PATH)

export_include: 
	cp -f *.h -t $(CONF_HEADERS_PATH)

clean:
	find . -regex '.*\.\(~\|o\|so\)' | xargs rm -rf
