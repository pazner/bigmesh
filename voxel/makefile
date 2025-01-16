-include ../Make.user

SOURCES = $(wildcard *.cpp)
OBJECTS = $(patsubst %.cpp,%.o,$(SOURCES))

LIB = ../libbigmesh.a

CPPFLAGS = -std=c++11 -O3 -g -I$(METIS_DIR)/include -I..
LDFLAGS = $(LIB) -L$(METIS_DIR)/lib -lmetis

APP = preprocess

.PHONY: clean

$(APP): $(OBJECTS)
	c++ $(CPPFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.cpp
	c++ $(CPPFLAGS) $< -c -o $@

clean:
	rm -f preprocess $(OBJECTS)
