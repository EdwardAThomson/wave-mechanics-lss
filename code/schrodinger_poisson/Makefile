CXX = g++
CXXFLAGS = -O2 -std=c++17 -Wall
LDFLAGS = -lfftw3 -lm

all: sp_1d sp_3d sp_3d_cosmo

sp_1d: sp_1d.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

sp_3d: sp_3d.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

sp_3d_cosmo: sp_3d_cosmo.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

clean:
	rm -f sp_1d sp_3d sp_3d_cosmo
	rm -rf output/

.PHONY: all clean
