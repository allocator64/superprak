CC = mpicxx
CFLAGS = -c -Wall -O2
LDFLAGS = 
SOURCES = main.cpp
OBJECTS = $(SOURCES:.cpp=.o)
FORMATTED = $(SOURCES:.cpp=.formatted)
BINARY = task2

all: $(FORMATTED) $(SOURCES) $(BINARY)

$(BINARY): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f $(OBJECTS) $(BINARY) $(FORMATTED)

.SUFFIXES: .formatted

.cpp.formatted:
	clang-format -style=Google $< >$@ && cp $@ $<
