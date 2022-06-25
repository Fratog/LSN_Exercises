CXXFLAGS := $(CXXFLAGS) -std=c++17 -O3 -Wall
CXX := g++

src_file := integral.cpp random.cpp
obj_file := $(patsubst %.cpp, %.o, $(src_file)) #alternatively $(src_file:.cpp=.o)
EXECUTABLE := integral.x

#$(patsubst %.o,%.c,$(objects)
##patsubst recognize the patterns % of file.o and modify it by adding $(OBJPATH) to the files given to the function //i believe so
VPATH=./SOURCE/

OBJPATH=./OBJECTS/


source := $(patsubst %.cpp, $(VPATH)%.cpp, $(src_file)) #i think here % works in another way, https://makefiletutorial.com/#-wildcard-1, in particular like a string subsittuion https://makefiletutorial.com/#string-substitution
objects := $(patsubst %.o, $(OBJPATH)%.o, $(obj_file)) 

#appo: 
#	@echo $(objects)	

all : $(objects)						#all : object file required, main.o random.o
	@echo $(obejcts)
	$(CXX) $(objects) -o $(EXECUTABLE) $(CXXFLAGS) #g++ .o files -o main.x flags
	
$(OBJPATH)%.o : %.cpp	
	#all file.o : the relative file.cpp https://makefiletutorial.com/#pattern-rules, here % works in pattern-rules mode, i guess
	#echo $(OBEJCTS)					
	$(CXX) -c $< -o $@ $(CXXFLAGS)  #$@--->targets--->$(OBJPATH)%.o  $<--first.cpp_in->%.cpp #g++ -c file.cpp -o relative file.o (do it for alle file.o)

.PHONY: clean
clean:
	rm $(OBJPATH)*.o $(EXECUTABLE) 
		 

##############
ifeq ("x","y")


mom: 
	echo $(src_file) # works
	echo $(obj_file) # works
	@echo "here we are " $(source) #//it seems it doesn't do anything
	@echo "here we are " $(objects)


# here's all your 'commented' makefile content...https://stackoverflow.com/questions/4493291/how-to-add-multi-line-comments-in-makefiles

appo : 
	echo $(OBJPATH)
	echo $(VPATH)

endif

#https://www.gnu.org/software/make/manual/html_node/Directory-Search.html#Directory-Search
