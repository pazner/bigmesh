# User must supply METIS_DIR in Make.user file
# Metis should be configured to use 64 bit integers
-include Make.user

SOURCES = $(wildcard *.cpp)
OBJECTS = $(patsubst %.cpp,%.o,$(SOURCES))

CPPFLAGS = -std=c++11 -O0 -g -I$(METIS_DIR)/include
LDFLAGS = -L$(METIS_DIR)/lib -lmetis

all: partition

partition: $(OBJECTS)
	c++ $(OBJECTS) $(LDFLAGS) -o partition

%.o: %.cpp
	c++ $(CPPFLAGS) $< -c -o $@

clean:
	rm -f $(OBJECTS)
	rm -f partition
