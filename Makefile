CXXFLAGS += -std=c++11 -g -Wall -Wextra
LDFLAGS = -static
EXTRA_CXXFLAGS = -DNDEBUG -O3
EXTRA_LDFLAGS =
INTERACTIONS = \
    jellium3.cpp \
    jellium4.cpp \
    ipl.cpp \
    lj.cpp \
    ljg.cpp \
    harddisk.cpp \

# override stuff here
-include features.mk

CXXFLAGS += $(EXTRA_CXXFLAGS)
LDFLAGS += $(EXTRA_LDFLAGS)

MAKEFILES = \
    Makefile \
    features.mk \

SOURCES = \
    $(INTERACTIONS) \
    storage.cpp \
    chainrunner.cpp \
    io.cpp \
    paircorrel.cpp \
    tools.cpp \

OBJECTS = $(SOURCES:%.cpp=%.o)

BINARIES = \
    postlhc \

# default target
all: ts.mk.hpp dep $(BINARIES)
	if [ -x ./push ]; then ./push; else true; fi

# depend on Makefiles
ts.mk.hpp: $(MAKEFILES)
	@touch $@

# create features.mk if absent
features.mk:
	@touch $@

# automagic dependencies
dep: $(OBJECTS) $(MAKEFILES)
	$(CXX) $(CXXFLAGS) -MM $(SOURCES) >.depend
-include .depend

postlhc: main.o $(OBJECTS)
	$(CXX) $(LDFLAGS) -o $@ $< $(OBJECTS)

clean:
	rm -f $(BINARIES) *.o .depend

.PHONY: all clean dep
