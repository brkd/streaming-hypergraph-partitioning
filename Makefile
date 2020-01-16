all:
	git add .
	git commit
	git push -u origin master
        g++ -o main main.cpp -std=c++17 -O3

