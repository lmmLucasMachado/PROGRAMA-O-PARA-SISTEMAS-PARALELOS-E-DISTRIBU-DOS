# obrigatorio a existencia das pastas obj, src, inc e bin 

CC := gcc

CFLAGS := -Wall

LFLAG := -fopenmp
  
SRCFILES := $(wildcard src/*.cpp)

all: $(SRCFILES:src/%.cpp=obj/%.o)
	$(CC) $(CFLAGS) obj/*.o $(LFLAG) -o bin/saida

obj/%.o : src/%.cpp
	$(CC) $(CFLAGS) -c $< -o $@

.PHONY: clean
clean:
	rm -rf obj/*
	rm -rf bin/*

run:
	bin/saida
