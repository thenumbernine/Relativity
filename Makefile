DIST_FILENAME=relativity
DIST_TYPE=app

include ../Common/Base.mk
include ../Tensor/Include.mk
include ../Profiler/Include.mk
include ../Parallel/Include.mk

.PHONY: test
test: default
	cd test; lua test.lua all

.PHONY: setbase
setbase: default
	cd test; lua test.lua all setbase

