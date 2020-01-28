CC := ifort
#local directroy with source files
SRC_DIR := src
#local directory to store build files
BLD_DIR := bld
BIN_DIR := bin

SRC_FILES := $(wildcard src/*.f90) # collect all source files
OBJ_FILES := $(patsubst src/%.f90, bld/%.o, $(SRC_FILES))  # define names for object files
MOD_FILES := $(patsubst src/%.f90, %.mod, $(SRC_FILES))

#$(info SRC_DIR="$(SRC_DIR)")
#$(info BLD_DIR="$(BLD_DIR)")

#$(info SRC_FILES="$(SRC_FILES)")
#$(info OBJ_FILES="$(OBJ_FILES)")

vpath %.f90 $(SRC_DIR)

define make-goal
$1/%.o: %.f90
	${CC} -c $$< -o $$@
	
# add dependencies
# can this be done automatically?
$1/pk_data_types.o: $1/pk_kinds.o $1/pk_enumeration.o $1/pk_utility.o
$1/pk_master.o: $1/pk_runner.o
$1/pk_runner.o: $1/pk_data_types.o $1/pk_read_input.o $1/pk_write_output.o $1/pk_time_control.o $1/pk_methods.o

endef

.PHONY: all checkdirs # these will not be treated as targets

all: checkdirs bin/point_kinetics_exe clean

bin/point_kinetics_exe: $(OBJ_FILES)
	${CC} $^ -o $@

checkdirs: $(BLD_DIR) $(BIN_DIR)

$(BLD_DIR):
	@mkdir -p $@

$(BIN_DIR):
	@mkdir -p $@

clean:
	@rm -rf $(BLD_DIR)
	@rm -f $(MOD_FILES)
	
$(eval $(call make-goal,$(BLD_DIR)))
