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
	@rm -rf bin/*||:
	@rm -rf obj/*||:

run:
	@mpirun --hostfile hostfile -ppn $(ppn) -n $$(($(node) * $(ppn))) bin/saida||:

help:
	@echo "Para compilar o programa basta"||:
	@echo "make"||:

	@echo "Para rodar o programa e necessario especificar a quantidade de nos (node) e programas por no (ppn)"||:
	@echo "Caso queira roda 12 programas em 3 nos diferentes basta rodar:"||:
	@echo "make run node=3 ppn=12"||:

	@echo "E para limpar os objetos e o executavel"||:
	@echo "make clean"||:


# HostList
# gateway chococino 164.41.127.81
# cm2 164.41.127.82
# cm3 164.41.127.83
# cm4164.41.127.84
