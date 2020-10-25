CC := mpicc

CFLAGS := 

LNFLAG := 

SRCFILES := $(wildcard src/*.c)

all: $(SRCFILES:src/%.c=obj/%.o)
	$(CC) $(CFLAGS) obj/*.o $(LNFLAG) -o bin/saida

obj/%.o : src/%.c
	$(CC) $(CFLAGS) -c $< -o $@

.PHONY: clean
clean:
	rm -rf obj/*
	rm -rf bin/*

run:
ifeq ($(node), 1)
	mpirun -host "cm2" -np $(np) bin/saida
endif
ifeq ($(node), 2)
	mpirun -host "cm2, cm3" -np $(np) bin/saida
endif
ifeq ($(node), 3)
	mpirun -host "cm2, cm3, cm4" -np $(np) bin/saida
else
	echo "A variavel node precisa estar entre 1 e 3"
endif


# HostList
# gateway chococino 164.41.127.81
# cm2 164.41.127.82
# cm3 164.41.127.83
# cm4164.41.127.84
