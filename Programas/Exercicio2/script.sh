# $1 = Quantidade de Bits
# $2 = Quantidade de processos por nó
# $3 = Quantidade de nós

mpicc -D MODULO_NUM_BITS=$1 -D NUM_PROCESS=$(($2 * $3)) -D NUM_NODES=$3 vow_with_hash.c -o main

echo "cm2:$2" > hostfile

if [ "$3" -gt 1 ]; then
    echo "cm3:$2" >> hostfile
fi

if [ "$3" -gt 2 ]; then
    echo "cm4:$2" >> hostfile
fi

mpirun --hostfile hostfile -n $(($2 * $3)) ./main > ./resultados/$1Bits_$2Processos_$3Nos
