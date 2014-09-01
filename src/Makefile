RELATIVE_PATH= ..
include $(RELATIVE_PATH)/common.inc

# -------------------------------
# -- Library paths and modules --
# -------------------------------
LIB= AcerDET

SUBTARGETS= compile_core \
	compile_configuration \
	compile_io \
	compile_analyse

OBJECTS= AcerDET.o

# ----------------------
# -- Parts of library --
# ----------------------

OBJS= $(shell find . -regex '.*\.o')

# ---------------------------
# -- Compile whole project --
# ---------------------------

all: $(SUBTARGETS) $(OBJECTS) libs

# ---------------------------
# -- Create Shared Library --
# ---------------------------

lib$(LIB).so: $(OBJS)
	$(CC) -shared $(OBJS) -o $@

lib$(LIB).a: $(OBJS)
	ar cr lib$(LIB).a $(OBJS)

libs: lib$(LIB).so lib$(LIB).a

# --------------------------
# -- Compile global files --
# --------------------------

%.o : %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

# -------------------------
# -- Library sub modules --
# -------------------------

compile_core:
	$(MAKE) -C ./core

compile_configuration:
	$(MAKE) -C ./conf

compile_io: 
	$(MAKE) -C ./io
	
compile_analyse:
	$(MAKE) -C ./analyse


clean_core:
	$(MAKE) -C ./core clean

clean_configuration:
	$(MAKE) -C ./conf clean

clean_io: 
	$(MAKE) -C ./io clean
	
clean_analyse:
	$(MAKE) -C ./analyse clean


# -------------------
# -- Clean objects --
# -------------------

clean: clean_core clean_configuration clean_io clean_analyse
	find . -regex '.*\.\(~\a|\|o\|so\)' | xargs rm -rf
	find . -regex '.*\.so\..*' | xargs rm -rf
