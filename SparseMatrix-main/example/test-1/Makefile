# by R.Lie, Nov 01, 2002

include /home/cyliu/AFEPack/Make.global_options

source = $(wildcard *.cpp)
object = $(patsubst %.cpp, %.o, $(source))
LDFLAGS += -L/home/cyliu/AFEPack/library/lib -lAFEPack
CXXFLAGS += -I /home/cyliu/SparseMatrix/SparseMatrix-main/include/
LIBS += -lsparse -L /home/cyliu/SparseMatrix/SparseMatrix-main/lib/

all : main

%.o : %.cpp
	$(CXX) -c -o $@ $< $(CXXFLAGS)

main : $(object)
	$(CXX) -o $@ $(object) $(LDFLAGS) $(LIBS)

clean :
	-rm -rf $(object)
	-rm -rf main
	-rm -f *.[nes]
	-rm -f *.dx

.PHONY : default clean
