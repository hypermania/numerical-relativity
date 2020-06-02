program_NAME := main
program_C_SRCS := $(wildcard *.c)
program_CXX_SRCS := $(wildcard *.cpp)
program_H_SRCS := $(wildcard *.h)
program_C_OBJS := ${program_C_SRCS:.c=.o}
program_CXX_OBJS := ${program_CXX_SRCS:.cpp=.o}
program_OBJS := $(program_C_OBJS) $(program_CXX_OBJS)
program_INCLUDE_DIRS :=
program_LIBRARY_DIRS :=
program_LIBRARIES :=

CPPFLAGS += $(foreach includedir,$(program_INCLUDE_DIRS),-I$(includedir))
LDFLAGS += $(foreach librarydir,$(program_LIBRARY_DIRS),-L$(librarydir)) 
LDFLAGS += $(foreach library,$(program_LIBRARIES),-l$(library)) 

CPPFLAGS += -std=c++1y

CPPFLAGS += -Wall

CPPFLAGS += -msse4.2 -mpopcnt -mbmi -mbmi2 -mavx -mavx2 -pthread

CPPFLAGS += -O3


.PHONY: all clean distclean

all: $(program_NAME)

$(program_OBJS): $(program_H_SRCS)

$(program_NAME): $(program_OBJS)
	$(LINK.cc) $(program_OBJS) -o $(program_NAME)

clean:
	$(RM) $(program_NAME)
	$(RM) $(program_OBJS)
	$(RM) $(wildcard *~)

distclean: clean

