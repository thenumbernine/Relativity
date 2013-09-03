TARGETDIR := /data/local/bin/
CPPFLAGS := -Wall -O0 -g

default: install

test: install
	cd test; lua test.lua all

setbase: install
	cd test; lua test.lua all setbase

DEPS := $(shell find src -name "*.h")

relativity: src/relativity.cpp $(DEPS)
	g++ $(CPPFLAGS) -std=c++0x src/relativity.cpp -o relativity

install: relativity
	cp relativity $(TARGETDIR) 

all: relativity

clean:
	-rm relativity

.PHONY: default all install test setbase plot clean



