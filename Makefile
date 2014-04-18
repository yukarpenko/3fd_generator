
ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)

CXX           = g++
F             = gfortran
CXXFLAGS      = -fPIC -g
LD            = g++
LDFLAGS       = -g
SOFLAGS       = -shared
FFLAGS        = -lgfortran -g

LIBS          = $(SYSLIBS)
GLIBS         = $(SYSLIBS)

OBJS        = main.o test.o readpt.o
              
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
