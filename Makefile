CXX = g++
CXXFLAGS = -std=c++14 -Wall

SRC_FILES = $(wildcard *.cpp)
OBJ_FILES = $(SRC_FILES:.cpp=.o)

all: rasterizer

rasterizer: $(OBJ_FILES)
	$(CXX) $(CXXFLAGS) -o $@ $^

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

clean:
	rm -f rasterizer $(OBJ_FILES)
