DEPDIR := .d
$(shell mkdir -p $(DEPDIR) >/dev/null)
DEPFLAGS = -MT $@ -MMD -MP -MF $(DEPDIR)/$*.Td
CC = mpicxx
CFLAGS = -c -Wall -O2 $(DEPFLAGS)
LDFLAGS = 
SOURCES = $(wildcard *.cpp)
HEADERS = $(wildcard *.h)
OBJECTS = $(SOURCES:.cpp=.o)
FORMATTED = $(SOURCES:.cpp=.cformatted) $(HEADERS:.h=.hformatted)
BINARY = task2
POSTCOMPILE = mv -f $(DEPDIR)/$*.Td $(DEPDIR)/$*.d

all: $(FORMATTED) $(SOURCES) $(BINARY)

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
