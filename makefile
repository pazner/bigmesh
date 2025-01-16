# User must supply METIS_DIR in Make.user file
# Metis should be configured to use 64 bit integers
-include Make.user

SOURCES = $(wildcard *.cpp)
OBJECTS = $(patsubst %.cpp,%.o,$(SOURCES))

CPPFLAGS = -std=c++11 -O3 -g -I$(METIS_DIR)/include
LDFLAGS = -L$(METIS_DIR)/lib -lmetis

LIB = libbigmesh.a

APPS = partition/partition

$(LIB): $(OBJECTS)
	ar -rcs $@ $^
	ranlib $(LIB)

$(APPS): %: %.cpp $(LIB)
	c++ $(CPPFLAGS) -I. -o $@ $< $(LIB) $(LDFLAGS)

%.o: %.cpp
	c++ $(CPPFLAGS) $< -c -o $@

clean:
	rm -f $(OBJECTS)
