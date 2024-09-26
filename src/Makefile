#------------------------------------------------------------------------
#
#  Some useful directories.
#
HOME   = $(shell cd ../../;pwd)
ABSDFTSRC=$(HOME)/dft/src
DFT    = ../
DFTSRC = ./
OML    = ../../oml

ifeq ($(OSTYPE),Linux)
	OPTFLAGS  =-O2 -fstrict-aliasing
	WARNFLAGS =-Wall -Winline 
  MAKE_COMMAND = make -j1
endif

ifeq ($(USER),hp)
	OPTFLAGS=-O2 -fomit-frame-pointer -fstrict-aliasing
	WARNFLAGS =-Wall -Winline
	MAKE_COMMAND = make -j2
endif

ifeq ($(OS),Windows_NT)
	OPTFLAGS  =-O3 -fstrict-aliasing -fno-keep-inline-functions
	WARNFLAGS =-Wall -Winline
	MAKE_COMMAND = make -j1
endif

#------------------------------------------------------------------------
#
#  Compile and link flags.
#
INCLUDE   =-I. -I$(OML)
DEBUG     =-O1 -DDEBUG 
NODEBUG   =-DNDEBUG $(OPTFLAGS)
CXXFLAGS  =$(INCLUDE) $(NODEBUG) $(WARNFLAGS) -fno-exceptions

#------------------------------------------------------------------------
#
#  What to compile
#
SOURCES := $(shell find . -name '*.C')
OBJECTS := $(SOURCES:.C=.o) $(wildcard repo/*.o)
DEPS    := $(SOURCES:.C=.d)

#------------------------------------------------------------------------
#
#  Now for the library targets
#
libdft.a : $(OBJECTS) 
	ar -rsc $@ $?

all : libdft.a omlstuff

OMLFLAGS='$(CXXFLAGS)'

omlstuff :
	-mkdir $(DFTSRC)/repo
	$(MAKE_COMMAND) -C $(OML) list TYPE='Polarization' PREFIX=Pol- REPO=$(ABSDFTSRC)/repo \
	DEFINES='-I$(ABSDFTSRC) -include $(ABSDFTSRC)/Misc/Polarization.H' \
	DEBUGFLAGS=$(OMLFLAGS)

%.d: %.C
	$(SHELL) -ec '$(CXX) -MM $(CXXFLAGS) $< | sed '\''s%$(*F).o%$(*D)/& $@%g'\'' > $@'

#	$(SHELL) -ec '$(CXX) -MM $(CXXFLAGS) $< | sed '/''s%$(*F).o%& $@%g'/'' > $@'
	

rmobjects :
	rm -f $(OBJECTS)
	
rmdep :
	rm -f $(DEPS)

rmjunk :
	find . -name '*~' -exec rm -f {} \;

checkdep :
	find . -name '*.d' -size 0 -print;
	find . -name '*.d' -size 0 -exec rm -f {} \;

#------------------------------------------------------------------------
#
#  Header dependancies
#

include $(DEPS)


