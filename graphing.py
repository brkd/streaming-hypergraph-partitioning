from os import listdir, stat
from os.path import isfile, join
from math import sqrt
from statistics import variance
import matplotlib.pyplot as plt
import threading
import subprocess
import sys
import csv

csv_header = ["Matrix", "rt(s)", "cuts", "slack", "imbal", "mb size", "hash", "part-count", "SD"]
def draw_graphs(partition_count, imbal, slack_val, matrix, mb_size, hash_count):
    X, Y_filler = [], []
    with open("N2Pvertex.txt") as n2pfile:
        for line in n2pfile:
            if '@' in line or '%' in line:
                continue
            else:        
                X.append(int(line.split(',')[0]))
                Y_filler.append(int(line.split(',')[2]))    
        Y = []
        for i in range(0, len(Y_filler)):
            if i == 0:
                Y.append(Y_filler[i])
            else:
                Y.append(Y[i - 1] + Y_filler[i])
    X2, Y2_filler = [], []
    with open("BFvertex.txt") as bffile:
        for line in bffile:
            if '@' in line or '%' in line:
                continue
            else:        
                X2.append(int(line.split(',')[0]))
                Y2_filler.append(int(line.split(',')[2]))    
        Y2 = []
        for i in range(0, len(Y2_filler)):
            if i == 0:
                Y2.append(Y2_filler[i])
            else:
                Y2.append(Y2[i - 1] + Y2_filler[i])
    X3, Y3_filler = [], []
    with open("BF2vertex.txt") as bf2file:
        for line in bf2file:
            if '@' in line or '%' in line:
                continue
            else:        
                X3.append(int(line.split(',')[0]))
                Y3_filler.append(int(line.split(',')[2]))    
        Y3 = []
        for i in range(0, len(Y3_filler)):
            if i == 0:
                Y3.append(Y3_filler[i])
            else:
                Y3.append(Y3[i - 1] + Y3_filler[i])
    X4, Y4_filler = [], []
    with open("BF3vertex.txt") as bf3file:
        for line in bf3file:
            if '@' in line or '%' in line:
                continue
            else:        
                X4.append(int(line.split(',')[0]))
                Y4_filler.append(int(line.split(',')[2]))    
        Y4 = []
        for i in range(0, len(Y4_filler)):
            if i == 0:
                Y4.append(Y3_filler[i])
            else:
                Y4.append(Y3[i - 1] + Y4_filler[i])
    plt.plot(X,Y,label='N2P')
    plt.plot(X2,Y2,label='BF1')
    plt.plot(X3,Y3,label='BF2')
    plt.plot(X4,Y4,label='BF3')
    plt.legend()
    plt.show()

def write_output(matrix_name, output, algorithm, partition_count, imbal, slack_val, randomization_count, byte_size = -1, hash_count = -1):
    if int(algorithm) == 2:
        name = "N2P"
    elif int(algorithm) == 3:
        name = "N2P_K"
    elif int(algorithm) == 4:
        name = "BF1"
    elif int(algorithm) == 5:
        name = "BF2"
    elif int(algorithm) == 6:
        name = "BF3"
    elif int(algorithm) == 7:
        name = "BF4MULTI"
    csv_name = "Results/" + name + ".csv"
    with open(csv_name, "a+") as f:
        writer = csv.writer(f)
        size_vecs = []
        size_vec = []
        durations = []
        cuts = []
        devs = []
        if stat(csv_name).st_size == 0:
            writer.writerow(csv_header)
        for line in output.splitlines():
            if "Duration:" in line:
                durations.append(line.split(":")[1])
            elif "Cuts:" in line:
                cuts.append(line.split(":")[1])
            elif "part " in line:
                size_vec.append(int(line.split(":")[1]))
                if len(size_vec) == int(partition_count):
                    size_vecs.append(size_vec)
                    size_vec = []
        for i in range(0, int(randomization_count)):
            devs.append(sqrt(variance(size_vecs[i])))
        for i in range(0, int(randomization_count)):
            writer.writerow([matrix_name, durations[i], cuts[i], slack_val, imbal, byte_size, hash_count, partition_count, devs[i], size_vecs[i]])
            
def start_bf_partitioning(algorithm, partition_count, imbal, slack_val, matrix, randomization_count, mb_size, hash_count):
    args = ("./main", algorithm, partition_count, imbal, slack_val, matrix, randomization_count, mb_size, hash_count)
    popen = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out = popen.stdout.read()
    if out != b'':
        write_output(matrix, out.decode('utf-8'), algorithm, partition_count, imbal, slack_val, randomization_count, mb_size, hash_count)

def start_partitioning(algorithm, partition_count, imbal, slack_val, matrix, randomization_count, k):
    if algorithm != 3:
        args = ("./main", algorithm, partition_count, imbal, slack_val, matrix, randomization_count)
    else:
        args = ("./main", algorithm, partition_count, imbal, slack_val, matrix, randomization_count, k)
    popen = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out = popen.stdout.read()
    if out != b'':
        write_output(matrix, out.decode('utf-8'), algorithm, partition_count, imbal, slack_val, randomization_count)

if __name__ == "__main__":
    if(len(sys.argv)) == 1:
        print("pc, imbal, slack_val, dir, mb_size, hc, n2p_k, rc")
        exit()
    elif sys.argv[1] == "-g" and len(sys.argv) == 2:
        print("pc, imbal, slack_val, matrix, mb_size, hc, n2p_k")
        exit()
    if sys.argv[1] == "-g" and len(sys.argv) > 8:
        algorithms = ["2", "3", "4", "5", "6"]
        partition_count = sys.argv[2]
        imbal = sys.argv[3]
        slack_val = sys.argv[4]
        matrix = sys.argv[5]
        mb_size = sys.argv[6]
        hash_count = sys.argv[7]
        k = sys.argv[8]
        threads = []
        for algorithm in algorithms:
            if int(algorithm) == 4 or int(algorithm) == 5 or int(algorithm) == 6:
                partitioner = threading.Thread(target=start_bf_partitioning, args=(algorithm, partition_count, imbal, slack_val, matrix, "1", mb_size, hash_count))
                threads.append(partitioner)
                partitioner.start()
            else:
                partitioner = threading.Thread(target=start_partitioning, args=(algorithm, partition_count, imbal, slack_val, matrix, "1", k))
                threads.append(partitioner)
                partitioner.start()
        for thread in threads:
            thread.join()
        draw_graphs(partition_count, imbal, slack_val, matrix, mb_size, hash_count)
    elif sys.argv[1] != "-g" and len(sys.argv) > 8:
        algorithms = ["2", "3", "4", "5", "6"]
        partition_count = sys.argv[1]
        imbal = sys.argv[2]
        slack_val = sys.argv[3]
        directory = sys.argv[4]
        mb_size = sys.argv[5]
        hash_count = sys.argv[6]
        k = sys.argv[7]
        rc = sys.argv[8]
        files = [f.split(".", 1)[0] for f in listdir(directory) if isfile(join(directory, f))]
        matrices = []
        for file in files:
            if file not in matrices:
                matrices.append(file)
                for matrix in matrices:
                    for algorithm in algorithms:
                        if int(algorithm) == 4 or int(algorithm) == 5 or int(algorithm) == 6:
                            partitioner = threading.Thread(target=start_bf_partitioning, args=(algorithm, partition_count, imbal, slack_val, directory + "/" + matrix, rc, mb_size, hash_count))
                            partitioner.start()
                        else:
                            partitioner = threading.Thread(target=start_partitioning, args=(algorithm, partition_count, imbal, slack_val, directory + "/" + matrix, rc, k))
                            partitioner.start()
                            
