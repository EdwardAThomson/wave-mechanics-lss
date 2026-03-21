CXX = g++
CXXFLAGS = -O2 -std=c++17 -Wall
LDFLAGS = -lfftw3 -lm

all: sp_1d

sp_1d: sp_1d.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

clean:
	rm -f sp_1d
	rm -rf output/

.PHONY: all clean
