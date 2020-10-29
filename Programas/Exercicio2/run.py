# nohup python3 run.py >& /dev/null &

import os
import time

def mais_de_24h(time):
    return time > 86400

bits = [24, 28, 31]

for b in bits:
    start = time.time()
    os.system('./script.sh {} {} {}'.format(b, 12, 3))
    end = time.time()
    if mais_de_24h(end-start):
        continue

    start = time.time()
    os.system('./script.sh {} {} {}'.format(b, 8, 3))
    end = time.time()
    if mais_de_24h(end-start):
        continue
    
    start = time.time()
    os.system('./script.sh {} {} {}'.format(b, 4, 3))
    end = time.time()
    if mais_de_24h(end-start):
        continue
    
    start = time.time()
    os.system('./script.sh {} {} {}'.format(b, 2, 3))
    end = time.time()
    if mais_de_24h(end-start):
        continue
    
    start = time.time()
    os.system('./script.sh {} {} {}'.format(b, 1, 3))
    end = time.time()
    if mais_de_24h(end-start):
        continue
    
    start = time.time()
    os.system('./script.sh {} {} {}'.format(b, 1, 2))
    end = time.time()
    if mais_de_24h(end-start):
        continue
    
    start = time.time()
    os.system('./script.sh {} {} {}'.format(b, 1, 1))
    end = time.time()
    if mais_de_24h(end-start):
        continue