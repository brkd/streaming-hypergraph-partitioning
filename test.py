import subprocess

args_1 = ("./main", "2", "16", "1.05", "100", "Matrices/wiki-Talk.mtx", "1")
args_2 = ("./main", "2", "2", "1.05", "100", "Matrices/crystk03.mtx", "1")
popen_1 = subprocess.Popen(args_1, stdout=subprocess.PIPE)
popen_2 = subprocess.Popen(args_2, stdout=subprocess.PIPE)

output_1 = popen_1.stdout.read()
output_2 = popen_2.stdout.read()
print(output_1)
print(output_2)
