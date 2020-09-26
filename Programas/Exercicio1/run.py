import os
threads = [1, 2, 4, 8, 12]
bits = [20, 24, 28]
runs = 10

for i in bits:
    for j in threads:
        os.system('./script.sh {} {} {}'.format(i, j, runs))
