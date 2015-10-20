
CXX           = g++
F             = gfortran
CXXFLAGS      = -fPIC -g -march=native $(shell root-config --cflags)
LD            = g++
LDFLAGS       = -g
FFLAGS        = -lgfortran -g -march=native

LIBS          = $(SYSLIBS) $(shell root-config --libs)
GLIBS         = $(SYSLIBS)

vpath %.cpp src
vpath %.f ../UKW
SRC        = main.cpp read.cpp generate.cpp cascade.cpp UKUtility.cpp \
     tree.cpp ParticlePDG2.cpp DecayChannel.cpp DatabasePDG2.cpp

SRC.UrQMD  = mainf.f upyth.f a.f siglookup.f addpart.f init.f output.f \
     string.f input.f tabinit.f paulibl.f angdis.f dectim.f \
     proppot.f ukw.f anndec.f delpart.f getmass.f upmerge.f detbal.f \
     getspin.f iso.f blockres.f ityp2pdg.f \
     dwidth.f jdecay2.f whichres.f boxprg.f hepchg.f \
     cascinit.f hepcmp.f make22.f hepnam.f  \
     coload.f numrec.f saveinfo.f scatter.f error.f

OBJS       = $(patsubst %.cpp,obj/%.o,$(SRC))

OBJS.UrQMD = $(patsubst %.f,obj/%.o,$(SRC.UrQMD))
              
TARGET	    = generator
#------------------------------------------------------------------------------
$(TARGET):       $(OBJS)  $(OBJS.UrQMD)
		$(LD) $(LDFLAGS) $^ -o $@ $(LIBS) $(FFLAGS)
		@echo "$@ done"
clean:
		@rm -f $(OBJS) $(OBJS.UrQMD) $(TARGET)

$(OBJS): | objdir

objdir:
	@mkdir -p obj

obj/%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@
	
obj/%.o : %.f
	$(F) $(FFLAGS) -c $< -o $@
