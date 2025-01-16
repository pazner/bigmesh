-include ../Make.user

SOURCES = $(wildcard *.cpp)
OBJECTS = $(patsubst %.cpp,%.o,$(SOURCES))

LIB = ../libbigmesh.a

CPPFLAGS = -std=c++11 -O0 -g -I$(METIS_DIR)/include -I..
LDFLAGS = -L$(METIS_DIR)/lib -lmetis $(LIB)

APP = preprocess

.PHONY: clean

$(APP): $(OBJECTS)
	c++ $(CPPFLAGS) $(LDFLAGS) -o $@ $^

%.o: %.cpp
	c++ $(CPPFLAGS) $< -c -o $@

clean:
	rm -f preprocess $(OBJECTS)
