SOURCES = $(wildcard *.cpp)
OBJECTS = $(patsubst %.cpp,%.o,$(SOURCES))

all: partition

partition: $(OBJECTS)
	c++ -std=c++11 $(OBJECTS) -o partition

%.o: %.cpp
	c++ -O0 -g -std=c++11 $< -c -o $@

clean:
	rm -f $(OBJECTS)
	rm -f partition
