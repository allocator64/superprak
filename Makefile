ifeq ($(OPENMP),1)
	CC = mpixlcxx_r -qsmp=omp -v
	USE_OPENMP = -D"USE_OPENMP"
else
	CC = mpicxx -v
	USE_OPENMP =
endif

CFLAGS = -O2 $(USE_OPENMP)
SOURCES = $(wildcard *.cpp)
HEADERS = $(wildcard *.h)
OBJECTS = $(SOURCES:.cpp=.o)
FORMATTED = $(SOURCES:.cpp=.cformatted) $(HEADERS:.h=.hformatted)
BINARY = task2

all: $(BINARY)

format: $(FORMATTED) $(SOURCES)

$(BINARY): $(SOURCES) $(HEADERS)
	$(CC) $(CFLAGS) $(SOURCES) -o $@

clean:
	rm -f $(OBJECTS) $(BINARY) $(FORMATTED) stderr-* stdout-*

.SUFFIXES: .cformatted .hformatted

.cpp.cformatted:
	clang-format -style=Google $< >$@ && cp $@ $<

.h.hformatted:
	clang-format -style=Google $< >$@ && cp $@ $<

graph.png: graph.gnuplot output.dat
	gnuplot <graph.gnuplot >graph.png
