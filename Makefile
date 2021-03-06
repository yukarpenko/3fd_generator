
CXX           = g++
F             = gfortran
CXXFLAGS      = -fPIC -O3 -march=native $(shell root-config --cflags)
LD            = g++
LDFLAGS       = -O3
FFLAGS        = -lgfortran -O3 -march=native

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
