# 
# Compile examples
#

SRC=$(wildcard *.c)
OUT=$(patsubst %.c,%.out,$(SRC))

NAMES=$(patsubst %.c,%,$(SRC))

.PHONY: test build $(NAMES)

build: $(OUT) 

%.out : %.c
	@echo "Compiling $<"
	@gcc -Wall -g -O0 $< -I../../target/cbindgen/  -L../../target/release -lproj4rs_c -o $@

clean:
	rm *.out

export LD_LIBRARY_PATH=$(shell realpath ../../target/release)

test: $(NAMES)

$(NAMES):
	@echo "==== $@ ===="
	@./$@.out

