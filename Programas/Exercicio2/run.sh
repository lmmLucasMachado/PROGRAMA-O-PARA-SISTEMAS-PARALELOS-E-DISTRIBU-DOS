#!/bin/bash

####################################
#                                  #
# Script para executar o código.   #
#                                  #
####################################

if [ "$#" -lt 1 ]; then
    echo "É necessário informar a quantidade de processos. Exemplo:"
    echo "./run.sh 4"
    exit
fi

mpirun -np $1 ./main
