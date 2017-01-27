TMPDIR:=.tmp
$(shell mkdir -p $(TMPDIR) >/dev/null)
ifeq ($(OPENMP),1)
	CC = mpixlcxx_r -qsmp=omp -v -O2 -D"USE_OPENMP"
	BINARY = task2_omp
else
	CC = mpicxx -Wall -O2
	BINARY = task2_no_omp
endif

SOURCES = $(wildcard src/*.cpp)
HEADERS = $(wildcard src/*.h)
OBJECTS = $(SOURCES:.cpp=.o)
FORMATTED = $(addprefix $(TMPDIR)/, $(notdir $(SOURCES)) $(notdir $(HEADERS)))

all: $(BINARY)

format: $(FORMATTED) $(SOURCES)

$(BINARY): $(SOURCES) $(HEADERS)
	$(CC) $(SOURCES) -o $@

clean:
	rm -f $(OBJECTS) $(BINARY) $(FORMATTED) stderr-* stdout-*

$(TMPDIR)/% : src/%
	clang-format -style=Google $< >$@ && cp $@ $<

graph.png: graph.gnuplot output.dat
	gnuplot <graph.gnuplot >graph.png
