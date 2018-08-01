

lib:
	g++ -O3 -Wall -shared -std=c++11 -fPIC steinhardt_binding.cpp steinhardt.cpp -o steinhardt.so `python -m pybind11 --includes` `python-config --cflags --ldflags --libs python`
new:
	g++ -std=c++11 steinhardt.cpp main.cpp -o steinhardt	