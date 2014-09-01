RELATIVE_PATH= ../..
include $(RELATIVE_PATH)/common.inc

OBJECTS= Configuration.o 

all: $(OBJECTS)

%.o : %.cpp
	$(CC) $(CFLAGS) -c $< -o $@ -lStd

clean:
	find . -regex '.*\.\(~\|o\|so\)' | xargs rm -rf
