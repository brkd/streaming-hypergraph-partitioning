all:
	echo algo pc imbal slack fn rc bs hc
	g++ -o main main.cpp -std=c++17 -O3 -fpermissive -fopenmp

