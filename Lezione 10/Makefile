CXXFLAGS := $(CXXFLAGS) -std=c++17 -O3 -Wall 
#CXXFLAGS := $(CXXFLAGS) -fno-elide-constructors -std=c++17 -O3 -Wall 

CXX := g++  
CXX2 := mpicxx
src_file :=  random.cpp chromosome.cpp population.cpp cities.cpp main.cpp #armadillo.cpp 
obj_file := $(patsubst %.cpp, %.o, $(src_file)) #alternatively $(src_file:.cpp=.o)

EXECUTABLE := main.out

VPATH=./SOURCE/
OBJPATH=./OBJECTS/

source := $(patsubst %.cpp, $(VPATH)%.cpp, $(src_file)) 
objects := $(patsubst %.o, $(OBJPATH)%.o, $(obj_file)) 


all : $(objects)						
	@echo $(obejcts)
	$(CXX2) $(objects) -o $(EXECUTABLE) $(CXXFLAGS) #g++ .o files -o main.out flags
	
$(OBJPATH)%.o : %.cpp	
	
	$(CXX2) -c $< -o $@ $(CXXFLAGS)  
	
launch : 
	mpiexec -np 8 $(EXECUTABLE)
	
move : 
	mv *.txt DATA/

.PHONY: clean 
clean:
	rm $(OBJPATH)*.o $(EXECUTABLE) 

###########################################

#https://www.gnu.org/software/make/manual/html_node/Directory-Search.html#Directory-Search
