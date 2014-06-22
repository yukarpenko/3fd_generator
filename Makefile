
CXX           = g++
F             = gfortran
CXXFLAGS      = -fPIC -g $(shell root-config --cflags)
LD            = g++
LDFLAGS       = -g $(shell root-config --libs)
FFLAGS        = -lgfortran -g

LIBS          = $(SYSLIBS)
GLIBS         = $(SYSLIBS)

OBJS        = main.o read.o ParticlePDG2.o DecayChannel.o DatabasePDG2.o
              
TARGET	    = generator
#------------------------------------------------------------------------------
$(TARGET):       $(OBJS)
		$(LD) $(LDFLAGS) $^ -o $@ $(LIBS) $(FFLAGS)
		@echo "$@ done"
clean:
		@rm -f $(OBJS) $(TARGET)

%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c $<
	
%.o : %.f
	$(F) $(FFLAGS) -c $<
