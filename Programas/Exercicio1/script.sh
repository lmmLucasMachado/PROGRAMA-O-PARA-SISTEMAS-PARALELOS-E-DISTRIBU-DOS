#!/bin/bash

#########################################################################
# Repete a execução do programa uma certa quantidade de vezes           #
# Ao começar a executar, pode sair do terminal ou fechar a conexão SSH  #
#                                                                       #                                                              #
# Para executar o script, use:                                          #
# nohup ./script.sh var1 var2 var3 >& /dev/null &                       #
#                                                                       #
# var1: nbits                                                           #
# var2: NUM_THREADS                                                     #
# var3: quantidade de execuções                                         #
#                                                                       #
# Não esqueça de modificar o programa com os valores corretos e         #
# compilar antes de executar o script                                   #
#                                                                       #
#########################################################################

if [ "$#" -lt 3 ]; then
    exit
fi

gcc -fopenmp -D MODULO_NUM_BITS=$1 -D NUM_THREADS=$2 -o main vow_with_hash.c

i=1

until [ $i -gt $3 ]
do
  echo -e "Run $i\n" >> resultados/$1bits$2threads.txt
  ./main >> resultados/$1bits$2threads.txt

  ((i=i+1))
done
