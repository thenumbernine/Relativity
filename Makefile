TARGETDIR := /data/local/bin/
CPPFLAGS := -Wall -O0

default: test

test: install
	cd test; lua test.lua

DEPS := $(shell find src -name "*.h")

relativity: src/relativity.cpp $(DEPS)
	g++ $(CPPFLAGS) -std=c++0x src/relativity.cpp -o relativity

install: relativity
	cp relativity $(TARGETDIR) 

all: relativity

clean:
	-rm relativity

.PHONY: default all install test plot clean



