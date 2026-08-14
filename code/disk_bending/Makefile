CXX      = g++
CXXFLAGS = -O3 -std=c++17 -Wall -Wextra -fopenmp
LDFLAGS  = -lfftw3 -lm -fopenmp

# Dependencies: FFTW3 only (LAPACK is deliberately not required; the
# tridiagonal eigensolver in src/eigen1d.h is self-contained and is certified
# against the analytic harmonic ladder by tests/test_harmonic.cpp).
#   Debian/Ubuntu: apt-get install libfftw3-dev
#   macOS/brew   : brew install fftw

HEADERS = $(wildcard src/*.h)

# The validation ladder. Fast enough to run routinely.
TESTS   = test_poisson test_free test_harmonic test_jeans test_sheet \
          test_poisson_xz

# Stage 1 measurements: minutes, not seconds, so kept out of `make check`.
MEASURE = test_bending

all: sheet_1d sheet_2d $(TESTS)

sheet_1d: sheet_1d.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

sheet_2d: sheet_2d.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

measure: $(MEASURE)

test_%: tests/test_%.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

# The validation ladder, in order. Each catches a different class of bug, so
# run them in sequence and stop at the first failure.
check: $(TESTS)
	@set -e; for t in $(TESTS); do \
		echo "-------------------------------------------------------------"; \
		./$$t; \
	done
	@echo "-------------------------------------------------------------"
	@echo "all validation tests passed"

clean:
	rm -f sheet_1d sheet_2d $(TESTS) $(MEASURE)

distclean: clean
	rm -rf output/

.PHONY: all check clean distclean
