main:	main.cpp	t-star.hpp
	g++ main.cpp -o main -std=c++11 -Wno-deprecated -I./lib/ -pthread -O2 -larmadillo

.PHONY: clean
clean:
	rm -f main
