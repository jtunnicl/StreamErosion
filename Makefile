CC=g++
CFLAGS=-Werror -std=c++11
CCFLAGS=-I.
LDFLAGS=
SOURCES=main.cpp topo.cpp flowrouting.cpp streamerosion.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=stream

all: $(SOURCES) $(EXECUTABLE)
    
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $(CCFLAGS) -c $< -o $@

.PHONY: clean
clean:
	rm -f $(EXECUTABLE) $(OBJECTS)
