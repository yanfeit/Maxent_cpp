CC=g++
#CFLAGS=-c -Wall -I ~/boost_1_55_0 -O3
CFLAGS=-c -Wall -I /opt/eigen/ -O3
LDFLAGS=
SOURCES=Maxent.cc
OBJECTS=$(SOURCES:.cc=.o)
EXECUTABLE=Maxent

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cc.o:
	$(CC) $(CFLAGS) $< -o $@
