####################################
#
# Use "make 127mer=1" to make 127mer version
# Use "make 63mer=1" to make 63mer version(default)
#
###################################

CC=             gcc
CPP=			g++
GCCVERSIONMAJOR := $(shell expr `$(CC) -dumpversion | cut -f1 -d.` \>= 4)
GCCVERSIONMINOR := $(shell expr `$(CC) -dumpversion | cut -f2 -d.` \>= 5)
ifdef debug
CFLAGS=          -O0 -g -fomit-frame-pointer  -DDEBUG
else
CFLAGS=          -O3 -fomit-frame-pointer  #-mcrc32 #-DDEBUG
endif
DFLAGS=         
OBJS=		 subhash.o kmerhash.o newhash.o arc.o attachPEinfo.o bubble.o bubbleStat.o check.o classifyEdge.o compactEdge.o \
		concatenateEdge.o connect.o contig.o cutTipPreGraph.o cutTip_graph.o \
		darray.o dfib.o dfibHeap.o fib.o fibHeap.o \
		hashFunction.o kmer.o lib.o loadGraph.o loadPath.o \
		loadPreGraph.o  localAsm.o main.o map.o mem_manager.o \
		node2edge.o orderContig.o output_contig.o output_pregraph.o \
		output_scaffold.o pregraph.o prlHashCtg.o prlHashReads.o prlRead2Ctg.o \
		prlRead2path.o  prlReadPaths.o  prlReadFillGap.o read2scaf.o readInterval.o stack.o\
		readseq1by1.o scaffold.o searchPath.o seq.o splitReps.o \
		cutTip_graph2.o linearEdge.o  read2edge.o iterate.o
#OBJS=prlReadPaths.o
PROG=           SOAPdenovo_LR-63mer
INCLUDES=	-Iinc
SUBDIRS=    sparsePregraph 
LIBPATH=	
LIBS=      	-pthread -lm -lrt -lbam -lz -L./inc
EXTRA_FLAGS=

BIT_ERR = 0
ifeq (,$(findstring $(shell uname -m), x86_64 ppc64 ia64))
BIT_ERR = 1
endif

ifdef 127mer
CFLAGS += -DMER127
ifdef debug
PROG = SOAPdenovo_LR-127mer_debug
else
PROG = SOAPdenovo_LR-127mer
endif
else
CFLAGS += -DMER63
ifdef debug
PROG = SOAPdenovo_LR-63mer_debug
else
PROG = SOAPdenovo_LR-63mer
endif
endif

LINUX = 0
ifneq (,$(findstring Linux,$(shell uname)))
LINUX = 1
EXTRA_FLAGS += -Wl,--hash-style=both
endif

ifneq (,$(findstring $(shell uname -m), x86_64))
CFLAGS += -m64
endif

ifneq (,$(findstring $(shell uname -m), ia64))
CFLAGS += 
endif

ifneq (,$(findstring $(shell uname -m), ppc64))
CFLAGS += -mpowerpc64
endif

.SUFFIXES:.c .o

.c.o:
		@printf "Compiling $<...                             \r"; \
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< || echo "Error in command: $(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $<"

all:	SOAPdenovo

.PHONY:all clean install

envTest:
		@test $(BIT_ERR) != 1 || sh -c 'echo "Fatal: 64bit CPU and Operating System required!";false;'
#		@test $(GCCVERSIONMAJOR) == 1 || sh -c 'echo "GCC version lower than 4.5.0";false;'
#		@test $(GCCVERSIONMINOR) == 1 || sh -c 'echo "GCC version lower than 4.5.0";false;'

ifdef 127mer
ifdef debug
SOAPdenovo:	envTest $(OBJS)
		@cd sparsePregraph;make 127mer=1 debug=1 clean all;cd ..;
		@$(CPP) sparsePregraph/*.o $(CFLAGS) -g -o $(PROG) $(OBJS) $(LIBPATH) $(LIBS) $(ENTRAFLAGS)
		@printf "Linking...\r"
		@printf "$(PROG) compilation done.\n";
else
SOAPdenovo:	envTest $(OBJS)
		@cd sparsePregraph;make 127mer=1  clean all;cd ..;
		@$(CPP) sparsePregraph/*.o $(CFLAGS) -o $(PROG) $(OBJS) $(LIBPATH) $(LIBS) $(ENTRAFLAGS)
		@printf "Linking...\r"
		@printf "$(PROG) compilation done.\n";
endif
else
ifdef debug
SOAPdenovo:	envTest $(OBJS)
		@cd sparsePregraph;make 63mer=1 debug=1 clean all;cd ..;
		@$(CPP) sparsePregraph/*.o $(CFLAGS) -g -o $(PROG) $(OBJS) $(LIBPATH) $(LIBS) $(ENTRAFLAGS)
		@printf "Linking...\r"
		@printf "$(PROG) compilation done.\n";
else
SOAPdenovo:	envTest $(OBJS)
		@cd sparsePregraph;make 63mer=1  clean all;cd ..;
		@$(CPP) sparsePregraph/*.o $(CFLAGS) -o $(PROG) $(OBJS) $(LIBPATH) $(LIBS) $(ENTRAFLAGS)
		@printf "Linking...\r"
		@printf "$(PROG) compilation done.\n";
endif	
endif

clean:
		@rm -fr gmon.out *.o a.out *.exe *.dSYM $(PROG) *~ *.a *.so.* *.so *.dylib
		@printf "$(PROG) cleaning done.\n";

install:
		@cp $(PROG) ../../bin/
		@printf "$(PROG) installed at ../../bin/$(PROG)\n"
