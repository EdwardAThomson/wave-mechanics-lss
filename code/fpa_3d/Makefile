CXX = g++
CXXFLAGS = -std=c++17 -O2 -Wall
LDFLAGS = -lfftw3 -lm

fpa_3d: fpa_3d.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

run: fpa_3d
	./fpa_3d

clean:
	rm -f fpa_3d fpa_3d_D*.dat

.PHONY: run clean
