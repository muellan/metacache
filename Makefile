REL_ARTIFACT  = metacache
DBG_ARTIFACT  = metacache_debug
TEST_ARTIFACT = run-tests

#COMPILER     = nvcc
COMPILER     = g++
MPI_COMPILER = mpicxx
DIALECT      = -std=c++11
WARNINGS     = -Wall -Wextra -Wpedantic

REL_FLAGS   = $(INCLUDES) $(MACROS) $(DIALECT) -O3 $(WARNINGS)
DBG_FLAGS   = $(INCLUDES) $(MACROS) $(DIALECT) -O0 -g $(WARNINGS)
TEST_FLAGS  = $(INCLUDES) $(MACROS) $(DIALECT) -O3 $(WARNINGS)

REL_LDFLAGS  = -pthread -s
DBG_LDFLAGS  = -pthread 
TEST_LDFLAGS = -pthread -s


#--------------------------------------------------------------------
HEADERS = \
		  src/alignment.h \
          src/args_handling.h \
          src/args_parser.h \
          src/bitmanip.h \
          src/candidates.h \
          src/chunk_allocator.h \
          src/cmdline_utility.h \
          src/config.h \
          src/dna_encoding.h \
          src/filesys_utility.h \
          src/hash_dna.h \
          src/hash_family.h \
          src/hash_int.h \
          src/hash_multimap.h \
          src/io_error.h \
          src/io_serialize.h \
          src/io_options.h \
          src/modes.h \
          src/parallel_task_queue.h \
          src/print_info.h \
          src/print_results.h \
          src/sequence_io.h \
          src/sequence_view.h \
          src/sketch_database.h \
          src/stat_confusion.h \
          src/stat_moments.h \
          src/statistics.h \
          src/taxonomy.h \
          src/taxonomy_io.h \
          src/timer.h \
          src/version.h

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
          src/print_info.cpp \
          src/print_results.cpp \
          src/sequence_io.cpp \
          src/taxonomy_io.cpp

TEST_HEADERS = \

TEST_SOURCES = \
          test\tests.cpp \
          test/hash_multimap_test.cpp


#--------------------------------------------------------------------
REL_DIR  = build_release
DBG_DIR  = build_debug
TEST_DIR = build_test

PLAIN_SRCS = $(notdir $(SOURCES))
PLAIN_OBJS = $(PLAIN_SRCS:%.cpp=%.o)

PLAIN_TEST_SRCS = $(notdir $(TEST_SOURCES))	
PLAIN_TEST_OBJS = $(PLAIN_TEST_SRCS:%.cpp=%.o)

#--------------------------------------------------------------------
REL_OBJS  = $(PLAIN_OBJS:%=$(REL_DIR)/%)
DBG_OBJS  = $(PLAIN_OBJS:%=$(DBG_DIR)/%)
TEST_OBJS = $(PLAIN_TEST_OBJS:%=$(TEST_DIR)/%)

REL_COMPILE  = $(COMPILER) $(REL_FLAGS) -c $< -o $@
DBG_COMPILE  = $(COMPILER) $(DBG_FLAGS) -c $< -o $@
TEST_COMPILE = $(COMPILER) $(TEST_FLAGS) -c $< -o $@



#--------------------------------------------------------------------
# main targets
#--------------------------------------------------------------------
.PHONY: all clean 
	
release: $(REL_DIR) $(REL_ARTIFACT)

debug: $(DBG_DIR) $(DBG_ARTIFACT)

test: $(TEST_DIR) $(TEST_ARTIFACT)

all: release debug test

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
	$(COMPILER) -o $(REL_ARTIFACT) $(REL_OBJS) $(REL_LDFLAGS)

$(REL_DIR)/main.o : src/main.cpp src/modes.h 
	$(REL_COMPILE)

$(REL_DIR)/print_info.o : src/print_info.cpp $(HEADERS)
	$(REL_COMPILE)

$(REL_DIR)/print_results.o : src/print_results.cpp $(HEADERS)
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

$(REL_DIR)/taxonomy_io.o : src/taxonomy_io.cpp $(HEADERS)
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

$(DBG_ARTIFACT): $(DBG_OBJS)
	$(COMPILER) -o $(DBG_ARTIFACT) $(DBG_OBJS) $(DBG_LDFLAGS)

$(DBG_DIR)/main.o : src/main.cpp src/modes.h 
	$(DBG_COMPILE)

$(DBG_DIR)/print_info.o : src/print_info.cpp $(HEADERS)
	$(DBG_COMPILE)

$(DBG_DIR)/print_results.o : src/print_results.cpp $(HEADERS)
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

$(DBG_DIR)/taxonomy_io.o : src/taxonomy_io.cpp $(HEADERS)
	$(DBG_COMPILE)

$(DBG_DIR)/sequence_io.o : src/sequence_io.cpp src/sequence_io.h src/io_error.h
	$(DBG_COMPILE)

$(DBG_DIR)/filesys_utility.o : src/filesys_utility.cpp src/filesys_utility.h
	$(DBG_COMPILE)

$(DBG_DIR)/cmdline_utility.o : src/cmdline_utility.cpp src/cmdline_utility.h
	$(DBG_COMPILE)

$(DBG_DIR)/args_handling.o : src/args_handling.cpp src/args_handling.h src/args_parser.h src/filesys_utility.h
	$(DBG_COMPILE)



#--------------------------------------------------------------------
# tests
#--------------------------------------------------------------------
$(TEST_DIR):
	mkdir $(TEST_DIR) 

$(TEST_ARTIFACT): $(TEST_OBJS)
	$(COMPILER) -o $(TEST_ARTIFACT) $(TEST_OBJS) $(TEST_LDFLAGS)

# test specific sources
$(TEST_DIR)/tests.o : test/tests.cpp $(TEST_HEADERS)
	$(TEST_COMPILE)
	
$(TEST_DIR)/hash_multimap_test.o : test/hash_multimap_test.cpp src/hash_multimap.h
	$(TEST_COMPILE)
	
