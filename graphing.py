from os import listdir, stat
from os.path import isfile, join
from math import sqrt
from statistics import variance
import matplotlib.pyplot as plt
import threading
import subprocess
import sys
import csv

csv_header = ["Algorithm", "rt(s)", "avg_cut", "max_cut", "min_cut", "slack", "imbal", "mb size", "hash", "part-count", "max_SD", "min_SD", "run_count", "spec_param"]
output_files = ["RANDOMvertex.txt", "N2Pvertex.txt", "N2P_Kvertex.txt", "BFvertex.txt", "BF2vertex.txt", "BF3vertex.txt", "BF4MULTIvertex.txt"]
def draw_graphs(partition_count, imbal, slack_val, matrix, mb_size, hash_count):    
    for file_name in output_files:
        with open(file_name) as n2pfile:
            X, Y_filler = [], []
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
            plt.plot(X,Y,label=file_name.split("vertex"))[0]
    plt.legend()
    plt.show()

def write_output(matrix_name, output, algorithm, partition_count, imbal, slack_val, randomization_count, byte_size = -1, hash_count = -1, spec_param = -1):
    name = matrix_name.split('/')[len(matrix_name.split('/')) - 1].replace(".*", "")
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
                durations.append(float(line.split(":")[1]))
            elif "Cuts:" in line:
                cuts.append(int(line.split(":")[1]))
            elif "part " in line:
                size_vec.append(int(line.split(":")[1]))
                if len(size_vec) == int(partition_count):
                    size_vecs.append(size_vec)
                    size_vec = []
        for i in range(0, int(randomization_count)):
            try:
                devs.append(sqrt(variance(size_vecs[i])))
            except:
                print(algorithm, name)
        max_cut = max(cuts)
        max_sd = ('%.2f'%max(devs))
        min_cut = min(cuts)
        min_sd = ('%.2f'%min(devs))
        avg_cut = ('%.2f'%(sum(cuts) / len(cuts)))
        avg_rt = ('%.2f'%(sum(durations) / len(cuts)))
        writer.writerow([algorithm, avg_rt, avg_cut, max_cut, min_cut, slack_val, imbal, byte_size, hash_count, partition_count, max_sd, min_sd, randomization_count, spec_param])            
            
def start_bf_partitioning(algorithm, partition_count, imbal, slack_val, matrix, randomization_count, mb_size, hash_count, num_layers):
    if algorithm != "7":
        args = ("./main", algorithm, partition_count, imbal, slack_val, matrix, randomization_count, mb_size, hash_count)
    else:
        args = ("./main", algorithm, partition_count, imbal, slack_val, matrix, randomization_count, mb_size, hash_count, num_layers)
    popen = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out = popen.stdout.read()
    if out != b'':
        if algorithm != "7":
            write_output(matrix, out.decode('utf-8'), algorithm, partition_count, imbal, slack_val, randomization_count, mb_size, hash_count)
        else:
            write_output(matrix, out.decode('utf-8'), algorithm, partition_count, imbal, slack_val, randomization_count, mb_size, hash_count, num_layers)

def start_partitioning(algorithm, partition_count, imbal, slack_val, matrix, randomization_count, k):
    if algorithm != "3":
        args = ("./main", algorithm, partition_count, imbal, slack_val, matrix, randomization_count)
    if algorithm == "0":
        args = ("./main", algorithm, partition_count, matrix, randomization_count)
    else:
        args = ("./main", algorithm, partition_count, imbal, slack_val, matrix, randomization_count, k)
    popen = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out = popen.stdout.read()
    if out != b'':
        if algorithm != "3":
            write_output(matrix, out.decode('utf-8'), algorithm, partition_count, imbal, slack_val, randomization_count)
        else:
            write_output(matrix, out.decode('utf-8'), algorithm, partition_count, imbal, slack_val, randomization_count, spec_param=k)
        
if __name__ == "__main__":
    if(len(sys.argv)) == 1:
        print("pc, imbal, slack_val, dir, mb_size, hc, n2p_k, num_layers, rc")
        exit()
    elif sys.argv[1] == "-g" and len(sys.argv) == 2:
        print("pc, imbal, slack_val, matrix, mb_size, hc, n2p_k, num_layers")
        exit()
    if sys.argv[1] == "-g" and len(sys.argv) > 8:
        algorithms = ["2", "3", "4", "5", "6", "7"]
        partition_count = sys.argv[2]
        imbal = sys.argv[3]
        slack_val = sys.argv[4]
        matrix = sys.argv[5]
        mb_size = sys.argv[6]
        hash_count = sys.argv[7]
        k = sys.argv[8]
        num_layers = sys.argv[9]
        threads = []
        for algorithm in algorithms:
            if int(algorithm) == 4 or int(algorithm) == 5 or int(algorithm) == 6 or int(algorithm) == 7:
                partitioner = threading.Thread(target=start_bf_partitioning, args=(algorithm, partition_count, imbal, slack_val, matrix, "1", mb_size, hash_count, num_layers))
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
        algorithms = ["2", "3", "4", "5", "6", "7"]
        partition_count = sys.argv[1]
        imbal = sys.argv[2]
        slack_val = sys.argv[3]
        directory = sys.argv[4]
        mb_size = sys.argv[5]
        hash_count = sys.argv[6]
        k = sys.argv[7]
        num_layers = sys.argv[8]
        rc = sys.argv[9]
        files = [f.split(".", 1)[0] for f in listdir(directory) if isfile(join(directory, f))]
        matrices = []
        for file in files:
            if file not in matrices:
                matrices.append(file)
        for matrix in matrices:
            for algorithm in algorithms:
                if int(algorithm) == 4 or int(algorithm) == 5 or int(algorithm) == 6 or int(algorithm) == 7:
                    partitioner = threading.Thread(target=start_bf_partitioning, args=(algorithm, partition_count, imbal, slack_val, directory + matrix, rc, mb_size, hash_count, num_layers))
                    partitioner.start()
                else:
                    partitioner = threading.Thread(target=start_partitioning, args=(algorithm, partition_count, imbal, slack_val, directory + matrix, rc, k))
                    partitioner.start()
                    
