#-----------------------------------------------------------------------------

PROGRAM = AliExternalTrackParam

HDRS =   AliVMisc.h  AliVParticle.h  AliVTrack.h  AliVVertex.h\
        TStatToolkit.h AliTMinuitToolkit.h AliPDG.h AliPID.h AliDrawStyle.h AliSysInfo.h AliXRDPROOFtoolkit.h TTreeStream.h AliLog.h \
        AliExternalTrackParam.h AliHelix.h 

OBJS = AliVMisc.o  AliVParticle.o  AliVTrack.o  AliVVertex.o \
        TStatToolkit.o AliTMinuitToolkit.o AliPDG.o AliPID.o AliDrawStyle.o AliSysInfo.o AliXRDPROOFtoolkit.o TTreeStream.o AliLog.o \
        AliExternalTrackParam.o  AliHelix.o \
        dict_$(PROGRAM).o

#-----------------------------------------------------------------------------
ObjSuf = o
SrcSuf = cxx
ExeSuf = exe
DllSuf = so
OutPutOpt = -o

ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --libs) -lMinuit
ROOTGLIBS = $(shell root-config --glibs)

# Directories
HEADRDIR = .

# debug option: -g
# no exceptions: -fno-exceptions

CXX = g++
CXXFLAGS =  -Wall -fPIC -g -std=c++11
LD = g++
LDFLAGS =  -Wall -g
SOFLAGS = -shared
ARCHFLAGS =

#
CXXFLAGS += $(ROOTCFLAGS) -I$(HEADRDIR) -I$(ROOTSYS)/include
LIBS = $(ROOTLIBS)
GLIBS = $(ROOTGLIBS)


%/%.o: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c "$<" -o "$@"

$(PROGRAM).so: $(OBJS)
	@echo "Linking $(PROGRAM).so ..."
	$(LD) $(ARCHFLAGS) $(LDFLAGS) $(SOFLAGS) $(OBJS) $(ROOTLIBS) -o$(PROGRAM).so
	@echo "linking done"

dict_$(PROGRAM).cpp: $(HDRS)
	@echo "Generating dictionary ..."
	rootcint -f dict_$(PROGRAM).cpp -c -p -I$(HEADRDIR) $(HDRS) Linkdef.h

clean:
	@rm -f $(OBJS) $(PROGRAM).so dict_*
 

#------------------------------------------------------------------------
