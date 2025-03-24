#Makefile for linux


EIGEN_PATH := /usr/local/include
ZLIB_INCLUDE := /usr/local/include
ZLIB_LIB := /usr/local/lib
DEBUG := 

CC = /opt/homebrew/Cellar/gcc/13.2.0/bin/gcc-13
CXX = /opt/homebrew/Cellar/gcc/13.2.0/bin/g++-13

ifdef DEBUG
CFLAGS = -g -O0 -Wall
CXXFLAGS = -g -O0 -Wall
else
CFLAGS = -O3 -Wall
CXXFLAGS = -O3 -Wall
endif

CPPFLAGS = 
LDFLAGS = -Llibs
LIBS =  -lm -lz -lomp

objs = $(patsubst %.cpp,%.o,$(wildcard src/*.cpp))

.PHONY: all
all: smr

smr: $(objs)
	$(CXX) $(CXXFLAGS) $(objs) $(LDFLAGS) $(LIBS) -o $@

smr_static: $(objs)
	$(CXX) $(CXXFLAGS) $(objs) $(LDFLAGS) -static $(LIBS) -o $@

$(objs): %.o: %.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

clean:
	@rm -rf src/*.o
	@rm -rf smr
	@rm -rf smr_static


install:
	@echo The binary is under building directory named smr.

