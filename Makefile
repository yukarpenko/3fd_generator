
CXX           = g++
F             = gfortran
CXXFLAGS      = -fPIC -g $(shell root-config --cflags)
LD            = g++
LDFLAGS       = -g $(shell root-config --libs)
FFLAGS        = -lgfortran -g

LIBS          = $(SYSLIBS)
GLIBS         = $(SYSLIBS)

vpath %.cpp src
SRC        = main.cpp read.cpp ParticlePDG2.cpp DecayChannel.cpp DatabasePDG2.cpp

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
