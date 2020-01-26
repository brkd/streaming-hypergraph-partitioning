from os import listdir
from os.path import isfile, join
import subprocess
import sys

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

matrix_path = "/home/hyperstream/partitioning/Matrices/"

matrix_names = [f for f in listdir(matrix_path) if isfile(join(matrix_path, f))]




#args_1 = ("./main", "2", "16", "1.05", "100", "Matrices/wiki-Talk.mtx", "1")
#args_2 = ("./main", "2", "2", "1.05", "100", "Matrices/crystk03.mtx", "1")
#popen_1 = subprocess.Popen(args_1, stdout=subprocess.PIPE)
#popen_2 = subprocess.Popen(args_2, stdout=subprocess.PIPE)

#output_1 = popen_1.stdout.read()
#output_2 = popen_2.stdout.read()


