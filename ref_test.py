from os import listdir, stat
from os.path import isfile, join
import threading
import subprocess
import csv
import sys

partition_counts = ["16"]
ref_size = ["1"]

csv_header = ["algo", "pc", "imbal", "slack", "rc", "total_reg", "total_r1", "total_r2", "total_r3", "ref_size(MB)", "dur_reg", "dur_r1", "dur_r2", "dur_r3"]
def write_output(matrix, out, alg, pc, imbal, slack, rc, ref_size):
    name = matrix.split('/')[len(matrix.split('/')) - 1].replace(".*", "")
    csv_name = name + ".csv"
    algo = "N2P"
    with open (csv_name, "a+") as f:
        writer = csv.writer(f)
        reg_t = []
        reg_c = []
        r1_t = []
        r1_c = []
        r2_t = []
        r2_c = []
        r3_t = []
        r3_c = []
        if stat(csv_name).st_size == 0:
            writer.writerow(csv_header)
        for line in out.splitlines():
            if "Duration:" in line:
                reg_t.append(float(line.split(":")[1]))
            elif "Cuts:" in line:
                reg_c.append(int(line.split(":")[1]))
            elif "Duration_R1:" in line:
                r1_t.append(float(line.split(":")[1]))
            elif "Cuts_R1:" in line:
                r1_c.append(int(line.split(":")[1]))
            elif "Duration_R2:" in line:
                r2_t.append(float(line.split(":")[1]))
            elif "Cuts_R2:" in line:
                r2_c.append(int(line.split(":")[1]))
            elif "Duration_R3:" in line:
                r3_t.append(float(line.split(":")[1]))
            elif "Cuts_R3:" in line:
                r3_c.append(int(line.split(":")[1]))
        reg_time = '%.2f'%(sum(reg_t) / len(reg_t))
        reg_cut = '%.2f'%(sum(reg_c) / len(reg_c))
        r1_time = '%.2f'%(sum(r1_t) / len(r1_t))
        r1_cut = '%.2f'%(sum(r1_c) / len(r1_c))
        r2_time = '%.2f'%(sum(r2_t) / len(r2_t))
        r2_cut = '%.2f'%(sum(r2_c) / len(r2_c))
        r3_time = '%.2f'%(sum(r3_t) / len(r3_t))
        r3_cut = '%.2f'%(sum(r3_c) / len(r3_c))
        writer.writerow(["N2P", pc, imbal, slack, rc, reg_cut, r1_cut, r2_cut, r3_cut, ref_size, reg_time, r1_time, r2_time, r3_time])       


def start_testing(algorithm, partition_count, imbal, slack, rc, matrix, ref_size):
    args = ("./main", algorithm, partition_count, imbal, slack, matrix, rc, ref_size)
    popen = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out = popen.stdout.read()
    if out != b'':
        write_output(matrix, out.decode('utf-8'), algorithm, partition_count, imbal, slack, rc, ref_size)


if __name__ == '__main__':
    if len(sys.argv) == 1:
        print("-------------")
        print("arg 1: imbal")
        print("arg 2: slack")
        print("arg 3: run count") 
        print("arg 4: directory -- All matrices in the directory given will be used.")
        print("-------------")
        exit()
    algorithm = "2"
    imbal = sys.argv[1]
    slack = sys.argv[2]
    rc = sys.argv[3] 
    directory = sys.argv[4]
    
    files = [f.split(".", 1)[0] for f in listdir(directory) if isfile(join(directory, f))]
    matrices = []
    for f in files:
        if f not in matrices:
            matrices.append(f)

    for matrix in matrices:
        for pc in partition_counts:
            for rs in ref_size:
                partitioner = threading.Thread(target=start_testing, args=(algorithm, pc, imbal, slack, rc, directory + matrix, rs))
                partitioner.start()
