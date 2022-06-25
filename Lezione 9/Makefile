CXXFLAGS := $(CXXFLAGS) -std=c++2a -O3 -Wall 
#CXXFLAGS := $(CXXFLAGS) -fno-elide-constructors -std=c++17 -O3 -Wall 
CXX := g++  

src_file :=  random.cpp statistics.cpp chromosome.cpp main.cpp population.cpp cities.cpp #armadillo.cpp 
obj_file := $(patsubst %.cpp, %.o, $(src_file)) #alternatively $(src_file:.cpp=.o)
EXECUTABLE := main.x

VPATH=./SOURCE/
OBJPATH=./OBJECTS/


source := $(patsubst %.cpp, $(VPATH)%.cpp, $(src_file)) 
objects := $(patsubst %.o, $(OBJPATH)%.o, $(obj_file)) 


all : $(objects)						
	@echo $(obejcts)
	$(CXX) $(objects) -o $(EXECUTABLE) $(CXXFLAGS) #g++ .o files -o main.x flags
	
$(OBJPATH)%.o : %.cpp	
	
	$(CXX) -c $< -o $@ $(CXXFLAGS)  
	
launch : 
	./$(EXECUTABLE)
	
move : 
	mv *.txt DATA/

.PHONY: clean 
clean:
	rm $(OBJPATH)*.o $(EXECUTABLE) 


###########################################

#https://www.gnu.org/software/make/manual/html_node/Directory-Search.html#Directory-Search
