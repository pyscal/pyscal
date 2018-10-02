

lib:
	g++ -O3 -Wall -shared -std=c++11 -fPIC steinhardt_binding.cpp steinhardt.cpp -o _steinhardt.so `python -m pybind11 --includes` `python-config --cflags --ldflags --libs python`
lib3:
	g++ -O3 -Wall -shared -std=c++11 -fPIC steinhardt_binding.cpp steinhardt.cpp -o _steinhardt.so `python3 -m pybind11 --includes` `python-config --cflags --ldflags --libs python3`

new:
	g++ -std=c++11 steinhardt.cpp main.cpp -o steinhardt	
