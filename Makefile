REL_ARTIFACT      = metacache
DBG_ARTIFACT      = metacache_debug
PRF_ARTIFACT      = metacache_prf

REL_CUDA_ARTIFACT = metacache_gpu
DBG_CUDA_ARTIFACT = metacache_gpu_debug
PRF_CUDA_ARTIFACT = metacache_gpu_prf

COMPILER      = $(CXX)
CUDA_COMPILER = nvcc
DIALECT       = -std=c++14
WARNINGS      = -Wall -Wextra -Wpedantic
NVCC_WARNINGS = -Xcompiler="-Wall -Wextra"
OPTIMIZATION  = -O3
#-march native -fomit-frame-pointer
# CUB = -I<path-to-cub>
NVCC_FLAGS    = $(CUB) -arch=sm_70 -lineinfo --expt-relaxed-constexpr --extended-lambda

REL_FLAGS        = $(INCLUDES) $(MACROS) $(DIALECT) $(OPTIMIZATION) $(WARNINGS)
DBG_FLAGS        = $(INCLUDES) $(MACROS) $(DIALECT) -O0 -g $(WARNINGS)
PRF_FLAGS        = $(INCLUDES) $(MACROS) $(DIALECT) $(OPTIMIZATION) -g $(WARNINGS)

REL_CUDA_FLAGS   = $(NVCC_FLAGS) $(INCLUDES) $(MACROS) $(DIALECT) $(OPTIMIZATION) $(NVCC_WARNINGS)
DBG_CUDA_FLAGS   = $(NVCC_FLAGS) $(INCLUDES) $(MACROS) $(DIALECT) -O0 -g $(NVCC_WARNINGS)
PRF_CUDA_FLAGS   = $(NVCC_FLAGS) $(INCLUDES) $(MACROS) $(DIALECT) $(OPTIMIZATION) -g $(NVCC_WARNINGS)

REL_LDFLAGS      = -pthread -s
DBG_LDFLAGS      = -pthread
PRF_LDFLAGS      = -pthread

REL_CUDA_LDFLAGS = $(NVCC_FLAGS) -Xcompiler="-pthread -s"
DBG_CUDA_LDFLAGS = $(NVCC_FLAGS) -Xcompiler="-pthread"
PRF_CUDA_LDFLAGS = $(NVCC_FLAGS) -Xcompiler="-pthread"


#--------------------------------------------------------------------
HEADERS = \
          src/alignment.h \
          src/batch_processing.h \
          src/bitmanip.h \
          src/candidate_generation.h \
          src/candidate_structs.h \
          src/chunk_allocator.h \
          src/classification.h \
          src/classification_statistics.h \
          src/cmdline_utility.h \
          src/config.h \
          src/database.h \
          src/dna_encoding.h \
          src/filesys_utility.h \
          src/gpu_hashmap.cuh \
          src/gpu_hashmap_operations.cuh \
          src/gpu_result_processing.cuh \
          src/hash_dna.h \
          src/hash_int.h \
          src/hash_multimap.h \
          src/host_hashmap.h \
          src/io_error.h \
          src/io_options.h \
          src/io_serialize.h \
          src/matches_per_target.h \
          src/mode_query.h \
          src/modes.h \
          src/options.h \
          src/printing.h \
          src/query_batch.cuh \
          src/querying.h \
          src/sequence_batch.cuh \
          src/sequence_io.h \
          src/sequence_view.h \
          src/stat_combined.cuh \
          src/stat_combined.h \
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
          src/database.cpp \
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
          src/taxonomy_io.cpp

CUDA_SOURCES = \
          src/gpu_hashmap.cu \
          src/query_batch.cu \
          src/sequence_batch.cu \
          src/stat_combined.cu


#--------------------------------------------------------------------
REL_DIR           = build_release
DBG_DIR           = build_debug
PRF_DIR           = build_profile

PLAIN_SRCS        = $(notdir $(SOURCES))
PLAIN_CUDA_SRCS   = $(notdir $(CUDA_SOURCES))
PLAIN_OBJS        = $(PLAIN_SRCS:%.cpp=%.o)
PLAIN_CUDA_OBJS   = $(PLAIN_CUDA_SRCS:%.cu=%.o)

#--------------------------------------------------------------------
REL_OBJS          = $(PLAIN_OBJS:%=$(REL_DIR)/%)
DBG_OBJS          = $(PLAIN_OBJS:%=$(DBG_DIR)/%)
PRF_OBJS          = $(PLAIN_OBJS:%=$(PRF_DIR)/%)

REL_CUDA_OBJS     = $(PLAIN_CUDA_OBJS:%=$(REL_DIR)/%)
DBG_CUDA_OBJS     = $(PLAIN_CUDA_OBJS:%=$(DBG_DIR)/%)
PRF_CUDA_OBJS     = $(PLAIN_CUDA_OBJS:%=$(PRF_DIR)/%)

REL_COMPILE       = $(COMPILER) $(REL_FLAGS) -c $< -o $@
DBG_COMPILE       = $(COMPILER) $(DBG_FLAGS) -c $< -o $@
PRF_COMPILE       = $(COMPILER) $(PRF_FLAGS) -c $< -o $@

REL_CUDA_COMPILE  = $(CUDA_COMPILER) $(REL_CUDA_FLAGS) -c $< -o $@
DBG_CUDA_COMPILE  = $(CUDA_COMPILER) $(DBG_CUDA_FLAGS) -c $< -o $@
PRF_CUDA_COMPILE  = $(CUDA_COMPILER) $(PRF_CUDA_FLAGS) -c $< -o $@


#--------------------------------------------------------------------
# main targets
#--------------------------------------------------------------------
.PHONY: all clean

release: $(REL_DIR) $(REL_ARTIFACT)

debug: $(DBG_DIR) $(DBG_ARTIFACT)

profile: $(PRF_DIR) $(PRF_ARTIFACT)

gpu_release: MACROS += -DGPU_MODE
gpu_release: $(REL_DIR) $(REL_CUDA_ARTIFACT)

gpu_debug: MACROS += -DGPU_MODE
gpu_debug: $(DBG_DIR) $(DBG_CUDA_ARTIFACT)

gpu_profile: MACROS += -DGPU_MODE
gpu_profile: $(PRF_DIR) $(PRF_CUDA_ARTIFACT)

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

$(REL_CUDA_ARTIFACT): $(REL_OBJS) $(REL_CUDA_OBJS)
	$(CUDA_COMPILER) -o $(REL_CUDA_ARTIFACT) $(REL_OBJS) $(REL_CUDA_OBJS) $(REL_CUDA_LDFLAGS)

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

$(REL_DIR)/database.o : src/database.cpp $(HEADERS)
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

$(REL_DIR)/gpu_hashmap.o : src/gpu_hashmap.cu $(HEADERS)
	$(REL_CUDA_COMPILE)

$(REL_DIR)/sequence_batch.o : src/sequence_batch.cu $(HEADERS)
	$(REL_CUDA_COMPILE)

$(REL_DIR)/stat_combined.o : src/stat_combined.cu $(HEADERS)
	$(REL_CUDA_COMPILE)

$(REL_DIR)/query_batch.o : src/query_batch.cu $(HEADERS)
	$(REL_CUDA_COMPILE)


#--------------------------------------------------------------------
# debug (out-of-place build)
#--------------------------------------------------------------------
$(DBG_DIR):
	mkdir $(DBG_DIR)

$(DBG_ARTIFACT): $(DBG_OBJS)
	$(COMPILER) -o $(DBG_ARTIFACT) $(DBG_OBJS) $(DBG_LDFLAGS)

$(DBG_CUDA_ARTIFACT): $(DBG_OBJS) $(DBG_CUDA_OBJS)
	$(CUDA_COMPILER) -o $(DBG_CUDA_ARTIFACT) $(DBG_OBJS) $(DBG_CUDA_OBJS) $(DBG_CUDA_LDFLAGS)

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

$(DBG_DIR)/database.o : src/database.cpp $(HEADERS)
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

$(DBG_DIR)/gpu_hashmap.o : src/gpu_hashmap.cu $(HEADERS)
	$(DBG_CUDA_COMPILE)

$(DBG_DIR)/sequence_batch.o : src/sequence_batch.cu $(HEADERS)
	$(DBG_CUDA_COMPILE)

$(DBG_DIR)/stat_combined.o : src/stat_combined.cu $(HEADERS)
	$(DBG_CUDA_COMPILE)

$(DBG_DIR)/query_batch.o : src/query_batch.cu $(HEADERS)
	$(DBG_CUDA_COMPILE)


#--------------------------------------------------------------------
# profile (out-of-place build)
#--------------------------------------------------------------------
$(PRF_DIR):
	mkdir $(PRF_DIR)

$(PRF_ARTIFACT): $(PRF_OBJS)
	$(COMPILER) -o $(PRF_ARTIFACT) $(PRF_OBJS) $(PRF_LDFLAGS)

$(PRF_CUDA_ARTIFACT): $(PRF_OBJS) $(PRF_CUDA_OBJS)
	$(CUDA_COMPILER) -o $(PRF_CUDA_ARTIFACT) $(PRF_OBJS) $(PRF_CUDA_OBJS) $(PRF_CUDA_LDFLAGS)

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

$(PRF_DIR)/database.o : src/database.cpp $(HEADERS)
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

$(PRF_DIR)/gpu_hashmap.o : src/gpu_hashmap.cu $(HEADERS)
	$(PRF_CUDA_COMPILE)

$(PRF_DIR)/sequence_batch.o : src/sequence_batch.cu $(HEADERS)
	$(PRF_CUDA_COMPILE)

$(PRF_DIR)/stat_combined.o : src/stat_combined.cu $(HEADERS)
	$(PRF_CUDA_COMPILE)

$(PRF_DIR)/query_batch.o : src/query_batch.cu $(HEADERS)
	$(PRF_CUDA_COMPILE)
