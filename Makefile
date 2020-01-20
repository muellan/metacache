REL_ARTIFACT  = metacache
DBG_ARTIFACT  = metacache_debug
PRF_ARTIFACT  = metacache_prf

COMPILER     = $(CXX)
DIALECT      = -std=c++14
WARNINGS     = -Wall -Wextra -Wpedantic
OPTIMIZATION = -O3
#-march native -fomit-frame-pointer

REL_FLAGS   = $(INCLUDES) $(MACROS) $(DIALECT) $(OPTIMIZATION) $(WARNINGS)
DBG_FLAGS   = $(INCLUDES) $(MACROS) $(DIALECT) -O0 -g $(WARNINGS)
PRF_FLAGS   = $(INCLUDES) $(MACROS) $(DIALECT) $(OPTIMIZATION) -g $(WARNINGS)

REL_LDFLAGS  = -pthread -s
DBG_LDFLAGS  = -pthread
PRF_LDFLAGS  = -pthread


#--------------------------------------------------------------------
HEADERS = \
          src/alignment.h \
          src/batch_processing.h \
          src/bitmanip.h \
          src/candidates.h \
          src/chunk_allocator.h \
          src/classification.h \
          src/classification_statistics.h \
          src/cmdline_utility.h \
          src/config.h \
          src/dna_encoding.h \
          src/filesys_utility.h \
          src/hash_dna.h \
          src/hash_int.h \
          src/hash_multimap.h \
          src/io_error.h \
          src/io_options.h \
          src/io_serialize.h \
          src/matches_per_target.h \
          src/modes.h \
          src/options.h \
          src/printing.h \
          src/querying.h \
          src/sequence_io.h \
          src/sequence_view.h \
          src/sketch_database.h \
          src/stat_confusion.h \
          src/stat_moments.h \
          src/string_utils.h \
          src/taxonomy.h \
          src/taxonomy_io.h \
          src/timer.h \
          src/version.h

SOURCES = \
          src/classification.cpp \
          src/cmdline_utility.cpp \
          src/filesys_utility.cpp \
          src/main.cpp \
          src/mode_build.cpp \
          src/mode_help.cpp \
          src/mode_info.cpp \
          src/mode_merge.cpp \
          src/mode_query.cpp \
          src/options.cpp \
          src/printing.cpp \
          src/sequence_io.cpp \
          src/sketch_database.cpp \
          src/taxonomy_io.cpp


#--------------------------------------------------------------------
REL_DIR  = build_release
DBG_DIR  = build_debug
PRF_DIR  = build_profile

PLAIN_SRCS = $(notdir $(SOURCES))
PLAIN_OBJS = $(PLAIN_SRCS:%.cpp=%.o)

#--------------------------------------------------------------------
REL_OBJS  = $(PLAIN_OBJS:%=$(REL_DIR)/%)
DBG_OBJS  = $(PLAIN_OBJS:%=$(DBG_DIR)/%)
PRF_OBJS  = $(PLAIN_OBJS:%=$(PRF_DIR)/%)

REL_COMPILE  = $(COMPILER) $(REL_FLAGS) -c $< -o $@
DBG_COMPILE  = $(COMPILER) $(DBG_FLAGS) -c $< -o $@
PRF_COMPILE  = $(COMPILER) $(PRF_FLAGS) -c $< -o $@



#--------------------------------------------------------------------
# main targets
#--------------------------------------------------------------------
.PHONY: all clean

release: $(REL_DIR) $(REL_ARTIFACT)

debug: $(DBG_DIR) $(DBG_ARTIFACT)

profile: $(PRF_DIR) $(PRF_ARTIFACT)

all: release debug profile test

clean :
	rm -rf build_*
	rm -f *.exe
	rm -f *.json
	rm -f $(REL_ARTIFACT)
	rm -f $(DBG_ARTIFACT)
	rm -f $(PRF_ARTIFACT)


#--------------------------------------------------------------------
# release (out-of-place build)
#--------------------------------------------------------------------
$(REL_DIR):
	mkdir $(REL_DIR)

$(REL_ARTIFACT): $(REL_OBJS)
	$(COMPILER) -o $(REL_ARTIFACT) $(REL_OBJS) $(REL_LDFLAGS)

$(REL_DIR)/main.o : src/main.cpp src/modes.h
	$(REL_COMPILE)

$(REL_DIR)/printing.o : src/printing.cpp $(HEADERS)
	$(REL_COMPILE)

$(REL_DIR)/classification.o : src/classification.cpp $(HEADERS)
	$(REL_COMPILE)

$(REL_DIR)/mode_build.o : src/mode_build.cpp $(HEADERS)
	$(REL_COMPILE)

$(REL_DIR)/mode_info.o : src/mode_info.cpp $(HEADERS)
	$(REL_COMPILE)

$(REL_DIR)/mode_merge.o : src/mode_merge.cpp $(HEADERS)
	$(REL_COMPILE)

$(REL_DIR)/mode_query.o : src/mode_query.cpp $(HEADERS)
	$(REL_COMPILE)

$(REL_DIR)/taxonomy_io.o : src/taxonomy_io.cpp $(HEADERS)
	$(REL_COMPILE)

$(REL_DIR)/sketch_database.o : src/sketch_database.cpp $(HEADERS)
	$(REL_COMPILE)

$(REL_DIR)/options.o : src/options.cpp $(HEADERS)
	$(REL_COMPILE)

$(REL_DIR)/mode_help.o : src/mode_help.cpp src/modes.h src/filesys_utility.h
	$(REL_COMPILE)

$(REL_DIR)/sequence_io.o : src/sequence_io.cpp src/sequence_io.h src/io_error.h
	$(REL_COMPILE)

$(REL_DIR)/filesys_utility.o : src/filesys_utility.cpp src/filesys_utility.h
	$(REL_COMPILE)

$(REL_DIR)/cmdline_utility.o : src/cmdline_utility.cpp src/cmdline_utility.h
	$(REL_COMPILE)


#--------------------------------------------------------------------
# debug (out-of-place build)
#--------------------------------------------------------------------
$(DBG_DIR):
	mkdir $(DBG_DIR)

$(DBG_ARTIFACT): $(DBG_OBJS)
	$(COMPILER) -o $(DBG_ARTIFACT) $(DBG_OBJS) $(DBG_LDFLAGS)

$(DBG_DIR)/main.o : src/main.cpp src/modes.h
	$(DBG_COMPILE)

$(DBG_DIR)/printing.o : src/printing.cpp $(HEADERS)
	$(DBG_COMPILE)

$(DBG_DIR)/classification.o : src/classification.cpp $(HEADERS)
	$(DBG_COMPILE)

$(DBG_DIR)/mode_build.o : src/mode_build.cpp $(HEADERS)
	$(DBG_COMPILE)

$(DBG_DIR)/mode_info.o : src/mode_info.cpp $(HEADERS)
	$(DBG_COMPILE)

$(DBG_DIR)/mode_merge.o : src/mode_merge.cpp $(HEADERS)
	$(DBG_COMPILE)

$(DBG_DIR)/mode_query.o : src/mode_query.cpp $(HEADERS)
	$(DBG_COMPILE)

$(DBG_DIR)/taxonomy_io.o : src/taxonomy_io.cpp $(HEADERS)
	$(DBG_COMPILE)

$(DBG_DIR)/sketch_database.o : src/sketch_database.cpp $(HEADERS)
	$(DBG_COMPILE)

$(DBG_DIR)/options.o : src/options.cpp $(HEADERS)
	$(DBG_COMPILE)

$(DBG_DIR)/mode_help.o : src/mode_help.cpp src/modes.h  src/filesys_utility.h
	$(DBG_COMPILE)

$(DBG_DIR)/sequence_io.o : src/sequence_io.cpp src/sequence_io.h src/io_error.h
	$(DBG_COMPILE)

$(DBG_DIR)/filesys_utility.o : src/filesys_utility.cpp src/filesys_utility.h
	$(DBG_COMPILE)

$(DBG_DIR)/cmdline_utility.o : src/cmdline_utility.cpp src/cmdline_utility.h
	$(DBG_COMPILE)


#--------------------------------------------------------------------
# profile (out-of-place build)
#--------------------------------------------------------------------
$(PRF_DIR):
	mkdir $(PRF_DIR)

$(PRF_ARTIFACT): $(PRF_OBJS)
	$(COMPILER) -o $(PRF_ARTIFACT) $(PRF_OBJS) $(PRF_LDFLAGS)

$(PRF_DIR)/main.o : src/main.cpp src/modes.h
	$(PRF_COMPILE)

$(PRF_DIR)/printing.o : src/printing.cpp $(HEADERS)
	$(PRF_COMPILE)

$(PRF_DIR)/classification.o : src/classification.cpp $(HEADERS)
	$(PRF_COMPILE)

$(PRF_DIR)/mode_build.o : src/mode_build.cpp $(HEADERS)
	$(PRF_COMPILE)

$(PRF_DIR)/mode_info.o : src/mode_info.cpp $(HEADERS)
	$(PRF_COMPILE)

$(PRF_DIR)/mode_merge.o : src/mode_merge.cpp $(HEADERS)
	$(PRF_COMPILE)

$(PRF_DIR)/mode_query.o : src/mode_query.cpp $(HEADERS)
	$(PRF_COMPILE)

$(PRF_DIR)/taxonomy_io.o : src/taxonomy_io.cpp $(HEADERS)
	$(PRF_COMPILE)

$(PRF_DIR)/sketch_database.o : src/sketch_database.cpp $(HEADERS)
	$(PRF_COMPILE)

$(PRF_DIR)/options.o : src/options.cpp $(HEADERS)
	$(PRF_COMPILE)

$(PRF_DIR)/mode_help.o : src/mode_help.cpp src/modes.h src/filesys_utility.h
	$(PRF_COMPILE)

$(PRF_DIR)/sequence_io.o : src/sequence_io.cpp src/sequence_io.h src/io_error.h
	$(PRF_COMPILE)

$(PRF_DIR)/filesys_utility.o : src/filesys_utility.cpp src/filesys_utility.h
	$(PRF_COMPILE)

$(PRF_DIR)/cmdline_utility.o : src/cmdline_utility.cpp src/cmdline_utility.h
	$(PRF_COMPILE)
