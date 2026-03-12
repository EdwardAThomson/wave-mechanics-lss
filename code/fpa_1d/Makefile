CXX = g++
CXXFLAGS = -std=c++17 -O2 -Wall
LDFLAGS = -lfftw3 -lm

fpa_1d: fpa_1d.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

run: fpa_1d
	./fpa_1d

clean:
	rm -f fpa_1d fpa_1d_D*.dat

.PHONY: run clean
