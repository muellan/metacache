REL_ARTIFACT = metacache
DBG_ARTIFACT = metacache_gdb

INCLUDES = 
#INCLUDES = -I stlsoft/include
MACROS   =

#COMPILER     = nvcc
COMPILER     = g++
MPI_COMPILER = mpicxx
DIALECT      = -std=c++11
WARNINGS     = -Wall -Wextra -Wpedantic

REL_FLAGS   = $(INCLUDES) $(MACROS) $(DIALECT) -O3 $(WARNINGS)
DBG_FLAGS   = $(INCLUDES) $(MACROS) $(DIALECT) -O0 -g $(WARNINGS)

REL_LDFLAGS =
DBG_LDFLAGS = 


#--------------------------------------------------------------------
HEADERS = src/alignment.h \
          src/args_parser.h \
          src/args_handling.h \
          src/bitmanip.h \
          src/cmdline_utility.h \
          src/dna_encoding.h \
          src/filesys_utility.h \
          src/hash_family.h \
          src/hash_multimap.h \
          src/io_error.h \
          src/min_hasher.h \
          src/modes.h \
          src/sequence_io.h \
          src/sequence_view.h \
          src/sketch_database.h \
          src/stat_moments.h \
          src/taxonomy.h \
          src/timer.h

SOURCES = \
    src/args_handling.cpp \
    src/cmdline_utility.cpp \
    src/filesys_utility.cpp \
    src/main.cpp \
    src/mode_annotate.cpp \
    src/mode_build.cpp \
    src/mode_help.cpp \
    src/mode_info.cpp \
    src/mode_query.cpp \
    src/sequence_io.cpp


#--------------------------------------------------------------------
REL_DIR = build_release
DBG_DIR = build_debug

PLAIN_SRCS = $(notdir $(SOURCES))	
PLAIN_OBJS = $(PLAIN_SRCS:%.cpp=%.o)

#--------------------------------------------------------------------
REL_OBJS = $(PLAIN_OBJS:%=$(REL_DIR)/%)
DBG_OBJS = $(PLAIN_OBJS:%=$(DBG_DIR)/%)

REL_COMPILE = $(COMPILER) $(REL_FLAGS) -c $< -o $@
DBG_COMPILE = $(COMPILER) $(DBG_FLAGS) -c $< -o $@

REL_LINK = $(COMPILER) -o $(REL_ARTIFACT) $(REL_OBJS) $(REL_LDFLAGS)
DBG_LINK = $(COMPILER) -o $(DBG_ARTIFACT) $(DBG_OBJS) $(DBG_LDFLAGS)



#--------------------------------------------------------------------
# main targets
#--------------------------------------------------------------------
.PHONY: all clean 
	
release: $(REL_DIR) $(REL_ARTIFACT)

debug: $(DBG_DIR) $(DBG_ARTIFACT)

all: release debug 

clean : 
	rm -rf build_*
	rm -f *.exe
	rm -f *.json
	rm -f $(REL_ARTIFACT)
	rm -f $(DBG_ARTIFACT)


#--------------------------------------------------------------------
# release
#--------------------------------------------------------------------
$(REL_DIR):
	mkdir $(REL_DIR) 

$(REL_ARTIFACT): $(REL_OBJS)
	$(REL_LINK)

$(REL_DIR)/main.o : src/main.cpp src/modes.h 
	$(REL_COMPILE)
	
$(REL_DIR)/mode_annotate.o : src/mode_annotate.cpp $(HEADERS)
	$(REL_COMPILE)
	
$(REL_DIR)/mode_build.o : src/mode_build.cpp $(HEADERS)
	$(REL_COMPILE)
	
$(REL_DIR)/mode_help.o : src/mode_help.cpp src/modes.h
	$(REL_COMPILE)

$(REL_DIR)/mode_info.o : src/mode_info.cpp $(HEADERS)
	$(REL_COMPILE)

$(REL_DIR)/mode_query.o : src/mode_query.cpp $(HEADERS)
	$(REL_COMPILE)

$(REL_DIR)/sequence_io.o : src/sequence_io.cpp src/sequence_io.h src/io_error.h
	$(REL_COMPILE)

$(REL_DIR)/filesys_utility.o : src/filesys_utility.cpp src/filesys_utility.h
	$(REL_COMPILE)

$(REL_DIR)/cmdline_utility.o : src/cmdline_utility.cpp src/cmdline_utility.h
	$(REL_COMPILE)

$(REL_DIR)/args_handling.o : src/args_handling.cpp src/args_handling.h src/args_parser.h src/filesys_utility.h
	$(REL_COMPILE)


#--------------------------------------------------------------------
# debug
#--------------------------------------------------------------------
$(DBG_DIR):
	mkdir $(DBG_DIR) 

$(DBG_DIR)/main.o : src/main.cpp src/modes.h
	$(DBG_COMPILE)
	
$(DBG_DIR)/mode_annotate.o : src/mode_annotate.cpp $(HEADERS)
	$(DBG_COMPILE)
	
$(DBG_DIR)/mode_build.o : src/mode_build.cpp $(HEADERS)
	$(DBG_COMPILE)
	
$(DBG_DIR)/mode_help.o : src/mode_help.cpp src/modes.h
	$(DBG_COMPILE)

$(DBG_DIR)/mode_info.o : src/mode_info.cpp $(HEADERS)
	$(DBG_COMPILE)

$(DBG_DIR)/mode_query.o : src/mode_query.cpp $(HEADERS)
	$(DBG_COMPILE)

$(DBG_DIR)/sequence_io.o : src/sequence_io.cpp src/sequence_io.h src/io_error.h
	$(DBG_COMPILE)

$(DBG_DIR)/filesys_utility.o : src/filesys_utility.cpp src/filesys_utility.h
	$(DBG_COMPILE)

$(DBG_DIR)/cmdline_utility.o : src/cmdline_utility.cpp src/cmdline_utility.h
	$(DBG_COMPILE)

$(DBG_DIR)/args_handling.o : src/args_handling.cpp src/args_handling.h src/args_parser.h src/filesys_utility.h
	$(DBG_COMPILE)
	
