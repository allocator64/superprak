DEPDIR := .d
$(shell mkdir -p $(DEPDIR) >/dev/null)
DEPFLAGS = -MT $@ -MMD -MP -MF $(DEPDIR)/$*.Td
ifeq ($(OPENMP),1)
	CC = mpicxx_r
	USE_OPENMP = "-DUSE_OPENMP"
else
	CC = mpicxx
	USE_OPENMP =
endif

CFLAGS = -c -Wall -O2 $(DEPFLAGS) $(USE_OPENMP)
LDFLAGS = 
SOURCES = $(wildcard *.cpp)
HEADERS = $(wildcard *.h)
OBJECTS = $(SOURCES:.cpp=.o)
FORMATTED = $(SOURCES:.cpp=.cformatted) $(HEADERS:.h=.hformatted)
BINARY = task2
POSTCOMPILE = mv -f $(DEPDIR)/$*.Td $(DEPDIR)/$*.d

all: $(BINARY)

format: $(FORMATTED) $(SOURCES)

$(BINARY): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

%.o : %.cpp
%.o : %.cpp $(DEPDIR)/%.d
	$(CC) $(CFLAGS) -o $@ $<
	$(POSTCOMPILE)


$(DEPDIR)/%.d: ;
.PRECIOUS: $(DEPDIR)/%.d

-include $(patsubst %,$(DEPDIR)/%.d,$(basename $(OBJECTS)))

clean:
	rm -f $(OBJECTS) $(BINARY) $(FORMATTED) stderr-* stdout-*

.SUFFIXES: .cformatted .hformatted

.cpp.cformatted:
	clang-format -style=Google $< >$@ && cp $@ $<

.h.hformatted:
	clang-format -style=Google $< >$@ && cp $@ $<

graph.png: graph.gnuplot output.dat
	gnuplot <graph.gnuplot >graph.png
