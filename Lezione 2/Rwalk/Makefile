CXXFLAGS := $(CXXFLAGS) -std=c++17 -O3 -Wall   
CXX := g++

src_file :=  random.cpp random_walk.cpp block_method.cpp  RW.cpp
obj_file := $(patsubst %.cpp, %.o, $(src_file)) 
EXECUTABLE := RW.x

VPATH=./SOURCE/
OBJPATH=./OBJECTS/

source := $(patsubst %.cpp, $(VPATH)%.cpp, $(src_file)) 
objects := $(patsubst %.o, $(OBJPATH)%.o, $(obj_file)) 

all : $(objects)						
	$(CXX) $(objects) -o $(EXECUTABLE) $(CXXFLAGS) 
	
$(OBJPATH)%.o : %.cpp	
	$(CXX) -c $< -o $@ $(CXXFLAGS)  

.PHONY: clean
clean:
	rm $(OBJPATH)*.o $(EXECUTABLE) 
		 

