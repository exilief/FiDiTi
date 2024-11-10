# Run: 'make' or 'make TARGET=Release'

TARGET := Debug
TARGETS := Debug Release

CXX ?= g++

CXXFLAGS_A := -std=c++17 -Wall
CXXFLAGS_D := $(CXXFLAGS_A) -g #-Og
CXXFLAGS_R := $(CXXFLAGS_A) -O3 -s -DNDEBUG
CPPFLAGS :=  # Preprocessor
LDFLAGS :=  # Linker

CXXFLAGS := $(CXXFLAGS_D)
Debug:   CXXFLAGS := $(CXXFLAGS_D)
Release: CXXFLAGS := $(CXXFLAGS_R)


# --------- #
# [Project] #
# --------- #

EXE := FiDiTi

SRC_DIRS := ./src ./include
BUILD_DIR := .
BIN_DIR := $(BUILD_DIR)/bin
OBJ_DIR := $(BUILD_DIR)/.obj



# ---------------- #
# [Compiler input] #
# ---------------- #

# Find all C / C++ source files
SRCS := $(shell find $(SRC_DIRS) -name '*.cpp' -or -name '*.c' -or -name '*.s')

# Prepend OBJ_DIR and append .o to every src file
OBJS = $(SRCS:%=$(OBJ_DIR)/%.o)

# For include dependencies (auto-generated)
DEPS := $(OBJS:.o=.d)

# Pass every folder in ./src to GCC so that it can find header files ('-I...')
INC_DIRS := $(shell find $(SRC_DIRS) -type d)
INC_FLAGS := $(addprefix -I,$(INC_DIRS))

# The -MMD and -MP flags generate Makefiles ('.d' extension, from cpp files)
CPPFLAGS += $(INC_FLAGS) -MMD -MP


# -------------- #
# [Main targets] #
# -------------- #

.PHONY: $(TARGET)

# Define target if TARGET is valid
Target := $(filter $(TARGET),$(TARGETS))

$(Target): $(BIN_DIR)/$(EXE)
	@echo Build: $(TARGET)

$(if $(Target),,$(TARGET)):
	@echo Error: invalid target '$(TARGET)'


# Main program
$(BIN_DIR)/$(EXE): $(OBJS)
	mkdir -p $(dir $@)
	$(CXX) $(OBJS) -o $@ $(LDFLAGS)


# Build step for C++ sources
$(OBJ_DIR)/%.cpp.o: %.cpp
	mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

# Build step for C sources
$(OBJ_DIR)/%.c.o: %.c
	mkdir -p $(dir $@)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@



clean:
	@echo Cleaning $(OBJ_DIR) $(BIN_DIR)
	rm -rf $(OBJ_DIR) $(BIN_DIR)

.PHONY: clean


# Include the .d makefiles. The - at the front suppresses the errors of missing Makefiles
# (Initially, all the .d files will be missing)
ifeq "" "$(strip $(filter clean,$(MAKECMDGOALS)))"
-include $(DEPS)
endif
