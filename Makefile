TARGETDIR := /data/local/bin/
CPPFLAGS := -Wall -O0

default: relativity

run: install
	cd tests; lua test.lua

DEPS := $(shell ls src/*.h)

relativity: src/relativity.cpp $(DEPS)
	g++ $(CPPFLAGS) -std=c++0x src/relativity.cpp -o relativity

install: relativity
	cp relativity $(TARGETDIR) 

all: relativity

clean:
	-rm relativity

.PHONY: default all install run plot clean



