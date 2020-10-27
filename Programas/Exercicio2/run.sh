#!/bin/bash

####################################
#                                  #
# Script para executar o código.   #
#                                  #
####################################

if [ "$#" -lt 1 ]; then
    echo "É necessário informar a quantidade de processos. Exemplo:"
    echo "./run.sh 4"
    echo "4 processos ao todo"
    echo "Não esqueça de mudar o hostfile com o número certo de processos por node, além de mudar o arquivo do programa"
    exit
fi

mpirun --hostfile hostfile -n $1 ./main