# ==============================
# Imbas_renderer Makefile
# Clean, modular, feature-driven
# ==============================

# -------- User options --------
# Optional external configuration (overrides defaults)
# Usage: make CONFIG=Config.mk or provide a custom file
CONFIG ?= Config.mk
-include $(CONFIG)

# Precedence: Makefile defaults < Config.mk < command line

SYSTEM ?= MacBook
ENABLE_OPENMP ?= 1
ENABLE_MPI    ?= 0

# -------- Compiler --------
CC ?= gcc-15

# -------- Target binary --------
BIN := imbas_renderer

# -------- Base flags --------
CFLAGS  := -Wall -Wextra -O3 -ffast-math -std=gnu17
LDFLAGS :=
LDLIBS  := -lm

# -------- Feature flags (OPT passed into code) --------

# -------- OpenMP --------
ifeq ($(ENABLE_OPENMP),1)
    OPT += -DENABLE_OPENMP

    ifeq ($(findstring clang,$(CC)),clang)
        CFLAGS  += -Xpreprocessor -fopenmp
        LDLIBS  += -lomp
    else
        CFLAGS  += -fopenmp
        LDFLAGS += -fopenmp
    endif
endif

# -------- MPI --------
ifeq ($(ENABLE_MPI),1)
    CC := mpicc
    OPT += -DENABLE_MPI
endif

# -------- Hybrid MPI + OpenMP notes --------
# If both ENABLE_MPI=1 and ENABLE_OPENMP=1:
# - Uses mpicc as compiler
# - Still applies OpenMP flags (-fopenmp or clang equivalent)
# This enables hybrid parallelism (MPI between nodes, OpenMP within node)

# -------- System-specific configs --------
ifeq ($(SYSTEM),MacBook)
    export MACOSX_DEPLOYMENT_TARGET := 26.4

    PNG_LIBS := /opt/homebrew/lib
    PNG_INCL := /opt/homebrew/include
    PNG_OPTS := -lpng

    HDF5_INCL := /opt/homebrew/include -DH5_USE_16_API
    HDF5_LIBS := /opt/homebrew/lib
    HDF5_OPTS := -lhdf5 -lz

    YAML_INCL := /opt/homebrew/include
    YAML_LIBS := /opt/homebrew/lib
    YAML_OPTS := -lyaml
endif

ifeq ($(SYSTEM),MacPro)
    PNG_LIBS := /usr/local/lib
    PNG_INCL := /usr/local/include
    PNG_OPTS := -lpng

    HDF5_INCL := /usr/local/include -DH5_USE_16_API
    HDF5_LIBS := /usr/local/lib
    HDF5_OPTS := -lhdf5 -lz

    YAML_INCL := /usr/local/include
    YAML_LIBS := /usr/local/lib
    YAML_OPTS := -lyaml
endif

ifeq ($(SYSTEM),Setonix)
    PNG_HOME := /software/projects/pawsey1164/cpower
    PNG_LIBS := $(PNG_HOME)/lib
    PNG_INCL := $(PNG_HOME)/include
    PNG_OPTS := -lpng
    HDF5_INCL := $(HDF5_DIR)/include -DH5_USE_16_API
    HDF5_LIBS := $(HDF5_DIR)/lib
    HDF5_OPTS := -lhdf5
endif

ifeq ($(SYSTEM),OzSTAR)
    PNG_LIBS := /opt/local/lib
    PNG_INCL := /opt/local/include
    PNG_OPTS := -lpng
    HDF5_INCL := $(EBROOTHDF5)/include -DH5_USE_16_API
    HDF5_LIBS := $(EBROOTHDF5)/lib
    HDF5_OPTS := -lhdf5
endif

ifeq ($(SYSTEM),Naranjo)
    PNG_HOME := /home/bsreedhar/software/modules/libpng/1.6.44
    PNG_LIBS := $(PNG_HOME)/lib
    PNG_INCL := $(PNG_HOME)/include
    PNG_OPTS := -lpng
    HDF5_INCL := $(EBROOTHDF5)/include -DH5_USE_16_API
    HDF5_LIBS := $(EBROOTHDF5)/lib
    HDF5_OPTS := -lhdf5
    YAML_INCL := 
    YAML_LIBS := 
    YAML_OPTS := -lyaml
endif

# ---------------------------------------------------------------------------
# Include paths
# ---------------------------------------------------------------------------
INCL := \
    -Isrc/core   \
    -Isrc/io     \
    -Isrc/tree   \
    -Isrc/sph    \
    -Isrc/render \
    -Isrc/parallel

ALL_INCL := $(INCL) \
    -I$(HDF5_INCL) \
    -I$(PNG_INCL)  \
    $(if $(YAML_INCL),-I$(YAML_INCL))

# ---------------------------------------------------------------------------
# Sources, objects, deps  (auto-discovered)
# ---------------------------------------------------------------------------
SRC_DIR   := src
BUILD_DIR := build

SRCS := $(shell find $(SRC_DIR) -name '*.c')
OBJS := $(patsubst $(SRC_DIR)/%.c, $(BUILD_DIR)/%.o, $(SRCS))
CFLAGS += $(OPT)
DEPS := $(OBJS:.o=.d)

# -------- Dependency include --------
# Include generated dependency files if they exist
-include $(DEPS)

# -------- Build rules --------
all: $(BIN)

$(BIN): $(OBJS)
	$(CC) -o $@ $^ \
	    -L$(HDF5_LIBS) $(HDF5_OPTS) \
	    -L$(PNG_LIBS)  $(PNG_OPTS)  \
	    $(if $(YAML_LIBS),-L$(YAML_LIBS)) $(YAML_OPTS) \
	    -lm

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(ALL_INCL) -MMD -MP -c $< -o $@

# -------- Utilities --------
clean:
	rm -rf $(BUILD_DIR) $(BIN)

print-config:
	@echo "SYSTEM        = $(SYSTEM)"
	@echo "CC            = $(CC)"
	@echo "ENABLE_OPENMP = $(ENABLE_OPENMP)"
	@echo "ENABLE_MPI    = $(ENABLE_MPI)"
	@echo "CFLAGS        = $(CFLAGS)"
	@echo "LDFLAGS       = $(LDFLAGS)"
	@echo "LDLIBS        = $(LDLIBS)"

.PHONY: all clean print-config
