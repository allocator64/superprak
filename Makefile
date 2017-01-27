TMPDIR:=.tmp
$(shell mkdir -p $(TMPDIR) >/dev/null)
ifeq ($(OPENMP),1)
	CC = mpixlcxx_r -qsmp=omp -v
	USE_OPENMP = -D"USE_OPENMP"
else
	CC = mpicxx
	USE_OPENMP =
endif

CFLAGS = -O2 $(USE_OPENMP)
SOURCES = $(wildcard src/*.cpp)
HEADERS = $(wildcard src/*.h)
OBJECTS = $(SOURCES:.cpp=.o)
FORMATTED = $(addprefix $(TMPDIR)/, $(notdir $(SOURCES)) $(notdir $(HEADERS)))
BINARY = task2

all: $(BINARY)

format: $(FORMATTED) $(SOURCES)

$(BINARY): $(SOURCES) $(HEADERS)
	$(CC) $(CFLAGS) $(SOURCES) -o $@

clean:
	rm -f $(OBJECTS) $(BINARY) $(FORMATTED) stderr-* stdout-*

$(TMPDIR)/%.cpp : src/%.cpp
	clang-format -style=Google $< >$@ && cp $@ $<

$(TMPDIR)/%.h : src/%.h
	clang-format -style=Google $< >$@ && cp $@ $<

graph.png: graph.gnuplot output.dat
	gnuplot <graph.gnuplot >graph.png
