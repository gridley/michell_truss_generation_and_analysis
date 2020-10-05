main: position.o
	g++ -g -std=c++17 -I/usr/include/eigen3 analyze_truss.cc position.o
position.o: position.cc position.h
	g++ -g -std=c++17 -c position.cc
