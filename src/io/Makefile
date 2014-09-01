RELATIVE_PATH= ../..
include $(RELATIVE_PATH)/common.inc

OBJECTS= InputRecord.o \
	OutputRecord.o

all: $(OBJECTS)

%.o : %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	find . -regex '.*\.\(~\|o\|so\)' | xargs rm -rf
