CXX = g++
CXXFLAGS = -std=c++17 -O2 -Wall
LDFLAGS = -lfftw3 -lm

all: fpa_3d_cosmo fpa_test_nu

fpa_3d_cosmo: src/fpa_3d_cosmo.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

fpa_test_nu: src/fpa_test_nu.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

run: fpa_3d_cosmo
	mkdir -p output/za
	./fpa_3d_cosmo

plot: run
	python3 scripts/plot_fpa_3d_cosmo.py

clean:
	rm -f fpa_3d_cosmo fpa_test_nu
	rm -rf output/*

.PHONY: all run plot clean
