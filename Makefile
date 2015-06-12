
CXX           = g++
F             = gfortran
CXXFLAGS      = -fPIC -g -march=native $(shell root-config --cflags) -fopenmp
LD            = g++
LDFLAGS       = -g -fopenmp
FFLAGS        = -lgfortran -g -march=native

LIBS          = $(SYSLIBS) $(shell root-config --libs)
GLIBS         = $(SYSLIBS)

vpath %.cpp src
SRC        = main.cpp read.cpp generate.cpp cascade.cpp UKUtility.cpp tree.cpp ParticlePDG2.cpp DecayChannel.cpp DatabasePDG2.cpp

OBJS       = $(patsubst %.cpp,obj/%.o,$(SRC)) 
              
TARGET	    = generator
#------------------------------------------------------------------------------
$(TARGET):       $(OBJS)
		$(LD) $(LDFLAGS) $^ -o $@ $(LIBS) $(FFLAGS)
		@echo "$@ done"
clean:
		@rm -f $(OBJS) $(TARGET)

$(OBJS): | objdir

objdir:
	@mkdir -p obj

obj/%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@
	
obj/%.o : %.f
	$(F) $(FFLAGS) -c $<
