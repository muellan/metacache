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
INCLUDE       =

NVCC_FLAGS    = -arch=$(CUDA_ARCH) -lineinfo --expt-relaxed-constexpr --extended-lambda
CXXFLAGS      = $(INCLUDE) $(MACROS) $(DIALECT) $(WARNINGS)

LDFLAGS       = -pthread

CUDA_FLAGS    = $(NVCC_FLAGS) $(INCLUDE) $(MACROS) $(DIALECT) $(NVCC_WARNINGS)
CUDA_LDFLAGS  = $(NVCC_FLAGS) -Xcompiler="-pthread"

# if MC_ZLIB=NO => deactivate zlib support
ifneq ($(MC_ZLIB),NO)
LDFLAGS += -lz
CUDA_LDFLAGS += -lz
else
MACROS += -DMC_NO_ZLIB
endif


#--------------------------------------------------------------------
HEADERS = \
          src/alignment.h \
          src/batch_processing.h \
          src/bitmanip.h \
          src/building.h \
          src/candidate_generation.h \
          src/candidate_structs.h \
          src/chunk_allocator.h \
          src/classification.h \
          src/classification_statistics.h \
          src/cmdline_utility.h \
          src/config.h \
          src/database.h \
          src/database_query.h \
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
          src/modes.h \
          src/options.h \
          src/printing.h \
          src/query_batch.cuh \
          src/query_handler.h \
          src/querying.h \
          src/sequence_batch.cuh \
          src/sequence_io.h \
          src/sequence_iostream.h \
          src/sequence_view.h \
          src/span.h \
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
          src/building.cpp \
          src/classification.cpp \
          src/cmdline_utility.cpp \
          src/database.cpp \
          src/filesys_utility.cpp \
          src/main.cpp \
          src/mode_build.cpp \
          src/mode_build_query.cpp \
          src/mode_help.cpp \
          src/mode_info.cpp \
          src/mode_merge.cpp \
          src/mode_query.cpp \
          src/options.cpp \
          src/printing.cpp \
          src/querying.cpp \
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

REL_CUDA_DIR      = build_gpu_release
DBG_CUDA_DIR      = build_gpu_debug
PRF_CUDA_DIR      = build_gpu_profile

PLAIN_SRCS        = $(notdir $(SOURCES))
PLAIN_CUDA_SRCS   = $(notdir $(CUDA_SOURCES))

PLAIN_OBJS        = $(PLAIN_SRCS:%.cpp=%.o)
PLAIN_CUDA_OBJS   = $(PLAIN_CUDA_SRCS:%.cu=%.o)

#--------------------------------------------------------------------
OBJS              = $(PLAIN_OBJS:%=$(DIR)/%)
CUDA_OBJS         = $(PLAIN_CUDA_OBJS:%=$(DIR)/%)

COMPILE           = $(COMPILER) $(CXXFLAGS) -c $< -o $@
CUDA_COMPILE      = $(CUDA_COMPILER) $(CUDA_FLAGS) -c $< -o $@


#--------------------------------------------------------------------
# main targets
#--------------------------------------------------------------------
release:
	$(MAKE) release_dummy DIR=$(REL_DIR) ARTIFACT=$(REL_ARTIFACT) MACROS=$(MACROS)

debug:
	$(MAKE) debug_dummy DIR=$(DBG_DIR) ARTIFACT=$(DBG_ARTIFACT) MACROS=$(MACROS)

profile:
	$(MAKE) profile_dummy DIR=$(PRF_DIR) ARTIFACT=$(PRF_ARTIFACT) MACROS=$(MACROS)

release_dummy: CXXFLAGS += $(OPTIMIZATION)
debug_dummy:   CXXFLAGS += -O0 -g
profile_dummy: CXXFLAGS += $(OPTIMIZATION) -g

release_dummy: LDFLAGS += -s

release_dummy: $(REL_DIR) $(REL_ARTIFACT)
debug_dummy:   $(DBG_DIR) $(DBG_ARTIFACT)
profile_dummy: $(PRF_DIR) $(PRF_ARTIFACT)

gpu_release:
	$(MAKE) gpu_release_dummy DIR=$(REL_CUDA_DIR) CUDA_ARTIFACT=$(REL_CUDA_ARTIFACT) MACROS="$(MACROS) -DGPU_MODE"

gpu_debug:
	$(MAKE) gpu_debug_dummy DIR=$(DBG_CUDA_DIR) CUDA_ARTIFACT=$(DBG_CUDA_ARTIFACT) MACROS="$(MACROS) -DGPU_MODE"

gpu_profile:
	$(MAKE) gpu_profile_dummy DIR=$(PRF_CUDA_DIR) CUDA_ARTIFACT=$(PRF_CUDA_ARTIFACT) MACROS="$(MACROS) -DGPU_MODE"

gpu_release_dummy: CXXFLAGS += $(OPTIMIZATION)
gpu_debug_dummy:   CXXFLAGS += -O0 -g
gpu_profile_dummy: CXXFLAGS += $(OPTIMIZATION) -g

gpu_release_dummy: LDFLAGS += -s

gpu_release_dummy: CUDA_FLAGS += $(OPTIMIZATION)
gpu_debug_dummy:   CUDA_FLAGS += -O0 -g
gpu_profile_dummy: CUDA_FLAGS += $(OPTIMIZATION) -g

gpu_release_dummy: CUDA_LDFLAGS += -Xcompiler="-s"

gpu_release_dummy: $(REL_CUDA_DIR) $(REL_CUDA_ARTIFACT)
gpu_debug_dummy:   $(DBG_CUDA_DIR) $(DBG_CUDA_ARTIFACT)
gpu_profile_dummy: $(PRF_CUDA_DIR) $(PRF_CUDA_ARTIFACT)


# phony targets
.PHONY: all clean gpu cpu
all: release debug profile
cpu: release
gpu: gpu_release

clean :
	rm -rf build_*
	rm -f *.exe
	rm -f *.json
	rm -f $(REL_ARTIFACT)
	rm -f $(DBG_ARTIFACT)
	rm -f $(PRF_ARTIFACT)
	rm -f $(REL_CUDA_ARTIFACT)
	rm -f $(DBG_CUDA_ARTIFACT)
	rm -f $(PRF_CUDA_ARTIFACT)


#--------------------------------------------------------------------
# out-of-place build
#--------------------------------------------------------------------
$(DIR):
	mkdir $(DIR)

$(ARTIFACT): $(OBJS)
	$(COMPILER) -o $(ARTIFACT) $(OBJS) $(LDFLAGS)

$(CUDA_ARTIFACT): $(OBJS) $(CUDA_OBJS)
	$(CUDA_COMPILER) -o $(CUDA_ARTIFACT) $(OBJS) $(CUDA_OBJS) $(CUDA_LDFLAGS)

$(DIR)/main.o : src/main.cpp src/modes.h
	$(COMPILE)

$(DIR)/building.o : src/building.cpp $(HEADERS)
	$(COMPILE)

$(DIR)/querying.o : src/querying.cpp $(HEADERS)
	$(COMPILE)

$(DIR)/printing.o : src/printing.cpp $(HEADERS)
	$(COMPILE)

$(DIR)/classification.o : src/classification.cpp $(HEADERS)
	$(COMPILE)

$(DIR)/mode_build.o : src/mode_build.cpp $(HEADERS)
	$(COMPILE)

$(DIR)/mode_build_query.o : src/mode_build_query.cpp $(HEADERS)
	$(COMPILE)

$(DIR)/mode_info.o : src/mode_info.cpp $(HEADERS)
	$(COMPILE)

$(DIR)/mode_merge.o : src/mode_merge.cpp $(HEADERS)
	$(COMPILE)

$(DIR)/mode_query.o : src/mode_query.cpp $(HEADERS)
	$(COMPILE)

$(DIR)/taxonomy_io.o : src/taxonomy_io.cpp $(HEADERS)
	$(COMPILE)

$(DIR)/database.o : src/database.cpp $(HEADERS)
	$(COMPILE)

$(DIR)/options.o : src/options.cpp $(HEADERS)
	$(COMPILE)

$(DIR)/mode_help.o : src/mode_help.cpp src/modes.h src/filesys_utility.h
	$(COMPILE)

$(DIR)/sequence_io.o : src/sequence_io.cpp src/sequence_io.h src/io_error.h src/sequence_iostream.h
	$(COMPILE)

$(DIR)/filesys_utility.o : src/filesys_utility.cpp src/filesys_utility.h
	$(COMPILE)

$(DIR)/cmdline_utility.o : src/cmdline_utility.cpp src/cmdline_utility.h
	$(COMPILE)

$(DIR)/gpu_hashmap.o : src/gpu_hashmap.cu $(HEADERS)
	$(CUDA_COMPILE)

$(DIR)/sequence_batch.o : src/sequence_batch.cu $(HEADERS)
	$(CUDA_COMPILE)

$(DIR)/stat_combined.o : src/stat_combined.cu $(HEADERS)
	$(CUDA_COMPILE)

$(DIR)/query_batch.o : src/query_batch.cu $(HEADERS)
	$(CUDA_COMPILE)
