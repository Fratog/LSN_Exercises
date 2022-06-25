CXXFLAGS := $(CXXFLAGS) -std=c++17 -O3 -Wall   
CXX := g++

src_file := options.cpp random.cpp block_method.cpp
obj_file := $(patsubst %.cpp, %.o, $(src_file)) 
EXECUTABLE := options.x

VPATH=./SOURCE/

OBJPATH=./OBJECTS/

source := $(patsubst %.cpp, $(VPATH)%.cpp, $(src_file)) 
objects := $(patsubst %.o, $(OBJPATH)%.o, $(obj_file)) 


all : $(objects)						
	$(CXX) $(objects) -o $(EXECUTABLE) $(CXXFLAGS) 
	
$(OBJPATH)%.o : %.cpp	
	$(CXX) -c $< -o $@ $(CXXFLAGS) 

launch : 
	./$(EXECUTABLE)
	
move : 
	mv *.txt DATA/

.PHONY: clean
clean:
	rm $(OBJPATH)*.o $(EXECUTABLE) 
		 


