#Makfile for linux

EIGEN_PATH := /usr/include
ZLIB_INCLUDE := /usr/include
ZLIB_LIB := /usr/lib64
DEBUG :=

CC = gcc
CXX = g++

ifdef DEBUG
CFLAGS = -g -O0 -Wall
CXXFLAGS = -g -O0 -Wall
else
CFLAGS = -O3 -Wall
CXXFLAGS = -O3 -Wall
endif

CPPFLAGS = -fopenmp
LDFLAGS = 
LIBS =  -lm -lz -lgomp

objs = $(patsubst %.cpp,%.o,$(wildcard src/*.cpp))

.PHONY: all
all: smr

smr: $(objs)
	$(CXX) $(CXXFLAGS) $(objs) $(LDFLAGS) $(LIBS) -o $@

smr_static: $(objs)
	$(CXX) $(CXXFLAGS) $(objs) $(LDFLAGS) -static $(LIBS) -pthread -ldl -o $@

$(objs): %.o: %.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

clean:
	@rm -rf src/*.o
	@rm -rf smr
	@rm -rf smr_static


install:
	@echo The binary is under building directory named smr.

