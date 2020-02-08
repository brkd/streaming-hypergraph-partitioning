from os import listdir
from os.path import isfile, join
import sys
import threading
import subprocess

imbal = 1.05
slack_val = 100
randomization_count = 1
byte_size = 10000000
hash_count = 4

def run_alg(matrix_name, alg, partition_count):
    if alg != 4:
        args = ("./main", alg, partition_count, imbal, slack_val, matrix_name, randomization_count)
        popen = subprocess.Popen(args, stdout=subprocess.PIPE)
    else:
        args = ("./main", alg, partition_count, imbal, slack_val, matrix_name, randomization_count, byte_size, hash_count)
        popen = subprocess.Popen(args, stdout=subprocess.PIPE)
    
def run_algs(dir, alg_array):
    partition_counts = [2, 16, 32, 1024, 4096]
    files = [f.split(".", 1)[0] for f in listdir(dir) if isfile(join(dir, f))]
    matrices = []

    for file in files:
        if file not in matrices:
            matrices.append(file)

     for alg in alg_array:
        for matrix in matrices:
            for pc in partition_counts:
                partitioning_thread = threading.Thread(target=run_alg, args=(matrix, alg, pc))
                partitioning_thread.start()

def run_all_algs(dir):
    alg_array = [1, 2, 3, 4]
    partition_counts = [2, 16, 32, 1024, 4096]
    files = [f.split(".", 1)[0] for f in listdir(dir) if isfile(join(dir, f))]
    matrices = []
    
    for file in files:
        if file not in matrices:
            matrices.append(file)
            
    for alg in alg_array:
        for matrix in matrices:
            for pc in partition_counts:
                partitioning_thread = threading.Thread(target=run_alg, args=(matrix, alg, pc))
                partitioning_thread.start()    
    

if __name__ == "__main__":
    dir = sys.argv[1]
    if len(sys.argv) == 2:
        run_all_algs(dir)
    else:
        run_algs(dir, sys.argv[2:])
            
        
        
