#-----------------------------------------------------------------------------

PROGRAM = AliSkimmedDataAnalysisMaker

HDRS =  TStatToolkit.h AliTMinuitToolkit.h AliPDG.h AliPID.h AliDrawStyle.h AliSysInfo.h AliXRDPROOFtoolkit.h TTreeStream.h AliLog.h AliExternalTrackParam.h V0FlatAna.h TrackFlatAna.h AliSkimmedDataAnalysisMaker.h\

OBJS = TStatToolkit.o AliTMinuitToolkit.o AliPDG.o AliPID.o AliDrawStyle.o AliSysInfo.o AliXRDPROOFtoolkit.o TTreeStream.o AliLog.o AliExternalTrackParam.o V0FlatAna.o TrackFlatAna.o AliSkimmedDataAnalysisMaker.o\
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
CXXFLAGS = -O3 -Wall -fPIC -g -std=c++11
LD = g++
LDFLAGS = -O3 -Wall -g 
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
