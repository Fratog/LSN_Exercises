CXXFLAGS := $(CXXFLAGS) -std=c++14 -O3 -Wall 
CXX := g++

src_file :=  random.cpp Monte_Carlo_ISING_1D.cpp
obj_file := $(patsubst %.cpp, %.o, $(src_file)) #alternatively $(src_file:.cpp=.o)
EXECUTABLE := Monte_Carlo_ISING_1D.x

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

.PHONY: clean clean_M clean_G
clean:
	rm $(OBJPATH)*.o $(EXECUTABLE) 
	#rm *.txt
	
clean_M:
	#rm DATA/Metropolis/*.dat
	#rm DATA/Eq_Metropolis/*.txt
	rm DATA/Metropolis/*.txt
	

clean_G:
#	rm DATA/Gibbs/*.dat
	#rm DATA/Eq_Gibbs/*.txt
	rm DATA/Gibbs/*.txt
	

