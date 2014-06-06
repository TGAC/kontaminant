#MAC=1
all		: kmer_contamination kmer_hash_build kmer_filter
#all	: cortex_con_cs solid_converter 
#all: solid_converter
ifndef CC
  CC = gcc	
endif

BIN = bin
LIB = lib

ifeq ($(MAXK),31)
   BITFIELDS = 1
endif

ifeq ($(MAXK),63)
   BITFIELDS = 2
endif

ifeq ($(MAXK),95)
   BITFIELDS = 3
endif

ifeq ($(MAXK),127)
   BITFIELDS = 4
endif

ifeq ($(MAXK),160)
   BITFIELDS = 5
endif

ifeq ($(MAXK),192)
   BITFIELDS = 6
endif

ifeq ($(MAXK),223)
   BITFIELDS = 7
endif

ifeq ($(MAXK),255)
   BITFIELDS = 8
endif

ifndef BITFIELDS
   BITFIELDS = 1
   MAXK = 31
endif 

# Main program includes
IDIR_BASIC =include/basic
IDIR_UTIL =include/util
IDIR_HASH  =include/hash_table 
IDIR_CORTEX_CON = include/cortex_con
IDIR_ALIGNMENT = include/alignment
IDIR_STATS = include/stats

# Test code includes
IDIR_BASIC_TESTS =include/test/basic
IDIR_HASH_TABLE_TESTS =include/test/hash_table
IDIR_CORTEX_CON_TESTS=include/test/cortex_con

#Default CUnit installation path in some Linux distributions. 
IDIR_CUNIT=/usr/local/include/CUnit/ 


UNAME := $(shell uname)


ifeq ($(UNAME), Darwin)
    MAC=1
endif

# Mac OS X specific. Assuming CUnit is installed with MacPorts
ifdef MAC
#¢-fnested-functions 
MACFLAG = -L/opt/local/lib/ 
IDIR_CUNIT=/opt/local/include/CUnit/ 
CFLAGS_CUNIT = -L/opt/local/lib/ -lncurses
#LD_CUNIT=
endif 

# 64bit architecture by default
ARCH =  -m64 
ifdef 32_BITS
 ARCH =  
endif

# Compiler options
OPT		= $(ARCH) -Wall -O3 $(MACFLAG) -DNUMBER_OF_BITFIELDS_IN_BINARY_KMER=$(BITFIELDS) -pthread -g

ifdef DEBUG
OPT	= $(ARCH) -Wall -O0 $(MACFLAG) -DNUMBER_OF_BITFIELDS_IN_BINARY_KMER=$(BITFIELDS) -g3 -pthread
endif

ifdef ENABLE_READ_PAIR
 OPT +=  -DENABLE_READ_PAIR
endif

ifdef ENABLE_READ_PAIR_OLD
 OPT += -DENABLE_READ_PAIR_OLD
endif

ifdef READ_PAIR_DEBUG_GRAPH
 OPT += -DREAD_PAIR_DEBUG_GRAPH
endif

ifdef DEBUG_PRINT_LABELS
 OPT += -DDEBUG_PRINT_LABELS
endif
kmerinfo:		OPT += -DNUMBER_OF_COLOURS=2
kmer_contamination:	OPT += -DKMER_CONTAMINATION -DNUMBER_OF_COLOURS=2 -DKMER_TOOLS
kmer_hash_build:	OPT += -DKMER_HASH_BUILD -DKMER_TOOLS	
kmer_filter:	OPT += -DKMER_TOOLS	
# Include dirs
CFLAGS_CORTEX_CON	= -I$(IDIR_BASIC) -I$(IDIR_HASH) -I$(IDIR_CORTEX_CON) -I$(IDIR_UTIL)
CFLAGS_CORTEX_BUB	= -I$(IDIR_BASIC) -I$(IDIR_HASH) -I$(IDIR_CORTEX_CON) -I$(IDIR_UTIL)
CFLAGS_GRAPHOUT         = -I$(IDIR_BASIC) -I$(IDIR_CORTEX_CON) -I$(IDIR_HASH) -I$(IDIR_UTIL)
CFLAGS_FILTERREADS         = -I$(IDIR_BASIC) -I$(IDIR_CORTEX_CON) -I$(IDIR_HASH) -I$(IDIR_UTIL)
CFLAGS_BUBBLEPARSE      = -I$(IDIR_BASIC) -I$(IDIR_CORTEX_CON) -I$(IDIR_HASH)
CFLAGS_KMERINFO         = -I$(IDIR_BASIC) -I$(IDIR_CORTEX_CON) -I$(IDIR_HASH)
CFLAGS_LOGPARSE         = -I$(IDIR_BASIC) -I$(IDIR_CORTEX_CON) -I$(IDIR_HASH) -I$(IDIR_UTIL)
CFLAGS_HASH_TABLE_TESTS	= -I$(IDIR_CUNIT) $(CFLAGS_CUNIT) -I$(IDIR_BASIC) -I$(IDIR_HASH) -I$(IDIR_CORTEX_CON) -I$(IDIR_HASH_TABLE_TESTS) 
CFLAGS_CORTEX_CON_TESTS	= -I$(IDIR_CUNIT) $(CFLAGS_CUNIT) -I$(IDIR_BASIC) -I$(IDIR_HASH) -I$(IDIR_CORTEX_CON) -I$(IDIR_CORTEX_CON_TESTS) 
CFLAGS_BASIC_TESTS	= -I$(IDIR_CUNIT) $(CFLAGS_CUNIT) -I$(IDIR_BASIC) -I$(IDIR_BASIC_TESTS)
CFLAGS_KMER_STATS   	= -I$(IDIR_BASIC) -I$(IDIR_CORTEX_CON) -I$(IDIR_HASH) -I$(IDIR_UTIL) -I$(IDIR_STATS)
CFLAGS_KMER_HASH_BUILD   	= -I$(IDIR_BASIC) -I$(IDIR_CORTEX_CON) -I$(IDIR_HASH) -I$(IDIR_UTIL) -I$(IDIR_STATS)
CFLAGS_DEMULTIPLEXER	= -I$(IDIR_BASIC) -I$(IDIR_UTIL) -I$(IDIR_STATS)
# Program objects
CORTEX_CON_OBJ = obj/kontaminant/file_format.o obj/kontaminant/flags.o obj/kontaminant/cleaning.o obj/kontaminant/path.o obj/kontaminant/perfect_path.o obj/kontaminant/branches.o obj/kontaminant/y_walk.o obj/kontaminant/cmd_line.o obj/kontaminant/binary_kmer.o obj/kontaminant/seq.o obj/kontaminant/element.o obj/kontaminant/hash_value.o obj/kontaminant/hash_table.o obj/kontaminant/dB_graph.o obj/kontaminant/file_reader.o obj/kontaminant/cortex_con.o obj/kontaminant/logger.o obj/kontaminant/metacortex.o obj/kontaminant/coverage_walk.o obj/util/node_queue.o

ifdef ENABLE_READ_PAIR
CORTEX_CON_OBJ += obj/kontaminant/binary_tree.o obj/kontaminant/read_pair.o
BIN_SUFFIX = $(join rp_,$(MAXK))
else
BIN_SUFFIX = $(MAXK)
endif 
cortex_con_rp:		CORTEX_CON_OBJ +=  obj/kontaminant/binary_tree.o obj/kontaminant/read_pair.o
cortex_con_mp:		CORTEX_CON_OBJ +=  obj/kontaminant/mark_pair.o
CORTEX_BUB_OBJ = $(CORTEX_CON_OBJ) obj/kontaminant/bubble_find.o

KMERINFO_OBJ =  obj/kontaminant/flags.o obj/kontaminant/path.o obj/kontaminant/binary_kmer.o obj/kontaminant/seq.o obj/kontaminant/element.o obj/kontaminant/hash_table.o obj/kontaminant/file_reader.o obj/util/kmerinfo.o obj/kontaminant/logger.o obj/kontaminant/hash_value.o
KMER_COUNT_OBJ = obj/kontaminant/flags.o  obj/kontaminant/binary_kmer.o  obj/kontaminant/path.o obj/kontaminant/perfect_path.o obj/kontaminant/seq.o obj/kontaminant/element.o obj/kontaminant/hash_value.o obj/kontaminant/hash_table.o obj/kontaminant/dB_graph.o obj/kontaminant/file_reader.o obj/stats/kmer_stats.o obj/stats/kmer_hash.o obj/kontaminant/logger.o obj/stats/kmer_reader.o obj/kontaminant/file_format.o
KMER_HASH_BUILD_OBJ = obj/kontaminant/flags.o obj/kontaminant/file_format.o obj/kontaminant/path.o obj/kontaminant/perfect_path.o obj/kontaminant/binary_kmer.o obj/kontaminant/seq.o obj/kontaminant/element.o obj/kontaminant/hash_value.o obj/kontaminant/hash_table.o obj/kontaminant/dB_graph.o obj/kontaminant/file_reader.o obj/stats/kmer_hash.o obj/stats/kmer_hash_build.o obj/stats/kmer_reader.o  obj/kontaminant/logger.o
KMER_FILTER_OBJ = obj/kontaminant/seq_io.o obj/kontaminant/flags.o obj/kontaminant/file_format.o  obj/kontaminant/path.o obj/kontaminant/perfect_path.o  obj/kontaminant/binary_kmer.o obj/kontaminant/seq.o obj/kontaminant/element.o obj/kontaminant/hash_value.o obj/kontaminant/hash_table.o obj/kontaminant/dB_graph.o obj/kontaminant/file_reader.o obj/stats/kmer_hash.o obj/stats/kmer_filter.o obj/stats/kmer_reader.o  obj/kontaminant/logger.o

#Library objects
LIBRARY_OBJ =  obj/kontaminant/file_format.o obj/kontaminant/analysis.o obj/kontaminant/flags.o obj/kontaminant/cleaning.o obj/kontaminant/path.o obj/kontaminant/perfect_path.o obj/kontaminant/branches.o obj/kontaminant/y_walk.o obj/kontaminant/cmd_line.o obj/kontaminant/binary_kmer.o obj/kontaminant/seq.o obj/kontaminant/element.o obj/kontaminant/hash_value.o obj/kontaminant/hash_table.o obj/kontaminant/dB_graph.o obj/kontaminant/file_reader.o obj/kontaminant/cortex_con.o obj/kontaminant/logger.o obj/kontaminant/metacortex.o obj/kontaminant/coverage_walk.o obj/util/node_queue.o

# Main rules


kmerinfo: remove_objects $(KMERINFO_OBJ)
	mkdir -p $(BIN); $(CC) -lm $(OPT) -o $(BIN)/kmerinfo $(KMERINFO_OBJ)


kmer_contamination: remove_objects $(KMER_COUNT_OBJ)
	mkdir -p $(BIN); $(CC) -lm $(OPT) -o $(BIN)/kmer_contamination_$(BIN_SUFFIX) $(KMER_COUNT_OBJ)

kmer_hash_build: remove_objects $(KMER_HASH_BUILD_OBJ)
	mkdir -p $(BIN); $(CC) -lm $(OPT) -o $(BIN)/kmer_hash_build_$(BIN_SUFFIX) $(KMER_HASH_BUILD_OBJ)

kmer_filter: remove_objects $(KMER_FILTER_OBJ)
	mkdir -p $(BIN); $(CC) -lm $(OPT) -o $(BIN)/kmer_filter_$(BIN_SUFFIX) $(KMER_FILTER_OBJ)



tests: remove_objects run_basic_tests run_hash_table_tests run_graph_tests

# Cleaning rules
.PHONY : clean
clean :
	rm -rf $(BIN)/*
	rm -rf obj

remove_objects:
	rm -rf obj

# Pattern rules

# cortex_con
obj/kontaminant/%.o : src/cortex_con/%.c include/cortex_con/%.h
	mkdir -p obj/kontaminant;  $(CC) $(CFLAGS_CORTEX_CON) $(OPT) -c $< -o $@

obj/kontaminant/%.o : src/basic/%.c include/basic/%.h
	mkdir -p obj/kontaminant;  $(CC) $(CFLAGS_CORTEX_CON) $(OPT) -c $< -o $@

obj/kontaminant/%.o : src/hash_table/hash_key/bob_jenkins/%.c include/hash_table/hash_value.h
	mkdir -p obj/kontaminant;  $(CC) $(CFLAGS_CORTEX_CON) $(OPT) -c $< -o $@

obj/kontaminant/%.o : src/hash_table/open_hash/%.c include/hash_table/open_hash/%.h
	mkdir -p obj/kontaminant;  $(CC) $(CFLAGS_CORTEX_CON) $(OPT) -c $< -o $@

obj/kontaminant/%.o : src/kontaminant/%.c 
	mkdir -p obj/kontaminant;  $(CC) $(CFLAGS_CORTEX_CON) $(OPT) -c $? -o $@

# Unit tests
obj/test/%.o : src/test/kontaminant/%.c include/test/kontaminant/%.h
	mkdir -p obj/test/; $(CC) $(CFLAGS_CORTEX_CON_TESTS) $(OPT) -c $< -o $@	

obj/test/run_dB_graph_tests.o : src/test/kontaminant/run_dB_graph_tests.c
	mkdir -p obj/test/; $(CC) $(CFLAGS_CORTEX_CON_TESTS) $(OPT) -c $< -o $@	

obj/test/%.o : src/basic/%.c include/basic/%.h
	mkdir -p obj/test/; $(CC) $(CFLAGS_BASIC_TESTS) $(OPT) -c $< -o $@

obj/test/%.o : src/test/basic/%.c include/test/basic/%.h
	mkdir -p obj/test/; $(CC) $(CFLAGS_BASIC_TESTS) $(OPT) -c $< -o $@

obj/test/%.o : src/test/hash_table/%.c include/test/hash_table/%.h
	mkdir -p obj/test/; $(CC) $(CFLAGS_HASH_TABLE_TESTS) $(OPT) -c $< -o $@

# Utils

obj/util/kmerinfo.o : src/util/kmerinfo.c 
	mkdir -p obj/util;  $(CC) $(CFLAGS_GRAPHOUT) $(OPT) -c $? -o $@


# Stats
obj/stats/%.o : src/stats/%.c include/stats/%.h
	mkdir -p obj/stats;  $(CC) $(CFLAGS_KMER_STATS) $(OPT) -c $< -o $@

obj/stats/kmer_stats.o : src/stats/kmer_stats.c
	mkdir -p obj/stats;  $(CC) $(CFLAGS_KMER_STATS) $(OPT) -c $< -o $@

obj/stats/demultiplexer.o : src/stats/demultiplexer.c
	mkdir -p obj/stats;  $(CC) $(CFLAGS_DEMULTIPLEXER) $(OPT) -c $< -o $@

obj/stats/kmer_hash_build.o : src/stats/kmer_hash_build.c
	mkdir -p obj/stats;  $(CC) $(CFLAGS_KMER_HASH_BUILD) $(OPT) -c $< -o $@
	
obj/stats/kmer_filter.o : src/stats/kmer_filter.c
	mkdir -p obj/stats;  $(CC) $(CFLAGS_KMER_HASH_BUILD) $(OPT) -c $< -o $@


libctx.dylib:$(LIBRARY_OBJ)
	ld  -dylib -dynamic $(LIBRARY_OBJ) -o $(LIB)/libctx_$(BIN_SUFFIX).dylib -lc -lz ; 

libctx.so.1:$(LIBRARY_OBJ)
	$(CC) -shared -Wl,-soname,libctx.so -o $(LIB)/libctx_$(BIN_SUFFIX).so.1 $(LIBRARY_OBJ) -lc -lz

dylib:
	mkdir -p $(BIN);
	@$(MAKE) clean; \
	case `uname` in \
		Linux) $(MAKE) CFLAGS="$(CFLAGS) -fPIC" libctx.so.1;; \
		Darwin) $(MAKE) CFLAGS="$(CFLAGS) -fPIC" libctx.dylib;; \
			*) echo 'Unknown OS';; \
	esac

