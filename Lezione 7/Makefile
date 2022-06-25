CXXFLAGS := $(CXXFLAGS) -std=c++14 -O3 -Wall  
CXX := g++

src_file := MD_MC.cpp random.cpp
obj_file := $(patsubst %.cpp, %.o, $(src_file)) 
EXECUTABLE := Molecular.x

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
	mv *.dat DATA/

.PHONY: clean clean_dat clean_frames clean_out
clean:
	rm $(OBJPATH)*.o $(EXECUTABLE) 

clean_dat:	
	rm  *.dat

clean_frames:	
	rm frames/solid/*.xyz	 

clean_out:
	rm *.out


