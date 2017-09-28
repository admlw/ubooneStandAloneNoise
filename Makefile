#makefile for gallery c++ programs.
#Note, being all-incllusive here: you can cut out libraries/includes you do not need
#you can also change the flags if you want too (Werror, pedantic, etc.)

CPPFLAGS=-I $(ROOT_INC) \

CXXFLAGS=-std=c++14 #-Wall #-pedantic #-Werror
CXX=g++

LDFLAGS=$$(root-config --libs) \

standalonenoise: standalonenoise.cc
		@$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $^

all: standalonenoise

clean:
	rm *.o standalonenoise

