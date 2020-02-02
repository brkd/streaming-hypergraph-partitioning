from os import listdir, stat
from os.path import isfile, join
from statistics import variance
from math import sqrt
import threading
import subprocess
import sys
import csv


csv_header = ["Algo", "Run-time", "Cuts", "Slack-val", "Imbal", "Byte-size", "Hash-count", "Part-count", "SD"]

def write_output(matrix_name, output, algorithm, partition_count, imbal, slack_value, randomization_count, byte_size, hash_count):
	csv_name = "Results/" + matrix_name.replace(".mtx", "") + ".csv"
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
			elif "part" in line:
				size_vec.append(int(line.split(":")[1]))
				if len(size_vec) == int(partition_count):
					size_vecs.append(size_vec)
					size_vec = []
		for i in range(0, int(randomization_count)):
			devs.append(sqrt(variance(size_vecs[i])))
		for i in range(0, int(randomization_count)):
			writer.writerow([algorithm, durations[i], cuts[i], slack_value, imbal, byte_size, hash_count, partition_count, devs[i], size_vecs[i]])

def start_partitioning(algorithm, partition_count, imbal, slack_value, matrix_name, randomization_count, byte_size=-1, hash_count=-1):
<<<<<<< HEAD
        args = ("./main", algorithm, partition_count, imbal, slack_value, "Matrices/" + matrix_name.replace(".mtx", ""), randomization_count)
        popen = subprocess.Popen(args, stdout=subprocess.PIPE)
        output = popen.stdout.read()
        write_output(matrix_name, output.decode('utf-8'), algorithm, partition_count, imbal, slack_value, randomization_count, byte_size, hash_count)
=======
	args = ("./main", algorithm, partition_count, imbal, slack_value, "/gandalf/data/Hyper/just_mtx/" + matrix_name, randomization_count)
	popen = subprocess.Popen(args, stdout=subprocess.PIPE)
	output = popen.stdout.read()
	write_output(matrix_name, output.decode('utf-8'), algorithm, partition_count, imbal, slack_value, randomization_count, byte_size, hash_count)
>>>>>>> 8e52ff8ee0dbe71cbb8b0ea3a4225e11b6de3e7e
	
if __name__ == "__main__":
	algorithm = sys.argv[1]
	partition_count = sys.argv[2]
	imbal = sys.argv[3]
	slack_value = sys.argv[4]
	randomization_count = sys.argv[5]
	if int(algorithm) == 4:
		if len(sys.argv) > 7:
			byte_size = sys.argv[6]
			hash_count = sys.argv[7]
		else:
			print("Info for BF partitioning not sufficient")

	matrix_path = "/home/brkd/Repos/streaming-hypergraph-partitioning/Matrices"

	matrix_names = [f for f in listdir(matrix_path) if isfile(join(matrix_path, f))]

	for matrix in matrix_names:
		if int(algorithm) == 4:
			partitioner = threading.Thread(target=start_partitioning, args=(algorithm, partition_count, imbal, slack_value, matrix, randomization_count, byte_size, hash_count))
			partitioner.start()
		else:
			partitioner = threading.Thread(target=start_partitioning, args=(algorithm, partition_count, imbal, slack_value, matrix, randomization_count))
			partitioner.start()

