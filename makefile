CXX       =g++
CXXFLAGS  =  -fpermissive -Wl,--no-as-needed -m64 -std=c++11
INCLUDE = -I include


#MKL------------------------------------------
MKL_TARGET    = mkl_test 
MKL_INCLUDE   =  -I ${MKLROOT}/include
MKL_INCLUDE += $(INCLUDE)

MKL_CXXFLAGS =
MKL_CXXFLAGS += $(CXXFLAGS)

MKL_LIBRARY   = -L ${MKLROOT}/lib/intel64

MKL_LDFLAGS  = -lmkl_rt -lpthread -lm -ldl 
#---------------------------------------------


#cuda------------------------------------------
CUDA_TARGET = cuda_test
CUDA_PATH ?= /usr/local/cuda-10.0
#CUDA_LDFLAGS = --dynamic-linker=/lib/ld-linux-armhf.so.3
#CUDA_LDFLAGS += $(addprefix -Xlinker ,$(CUDA_LDFLAGS))
#CUDA_LDFLAGS += --sysroot=$(TARGET_FS)
#CUDA_LDFLAGS += -rpath-link=$(TARGET_FS)/lib -L $(TARGET_FS)/lib
#CUDA_LDFLAGS += -rpath-link=$(TARGET_FS)/usr/lib -L $(TARGET_FS)/usr/lib
#CUDA_LDFLAGS += -rpath-link=$(TARGET_FS)/usr/lib/aarch64-linux-gnu -L $(TARGET_FS)/usr/lib/aarch64-linux-gnu
#CUDA_LDFLAGS += --unresolved-symbols=ignore-in-shared-libs


CUDA_CXXFLAGS = 
CUDA_CXXFLAGS += -std=c++11
CUDA_CXXFLAGS += -m64 
CUDA_CXXFLAGS += -g -G 
#CUDA_CXXFLAGS += -isystem=$(TARGET_FS)/usr/include
#CUDA_CXXFLAGS += -isystem=$(TARGET_FS)/usr/include/aarch64-linux-gnu
#CUDA_CXXFLAGS += $(addprefix -Xcompiler ,$(CUDA_CXXFLAGS))

#CUDA_CXXFLAGS += -mfloat-abi=hard
#CUDA_CXXFLAGS += $(CXXFLAGS)


NVCC          := $(CUDA_PATH)/bin/nvcc -ccbin $(CXX)
CUDA_INCLUDE  = -I include/cuda_inc
CUDA_INCLUDE += $(INCLUDE) 

CUDA_LDFLAGS = -lcublas
#----------------------------------------------

OBJ_DIR = ./obj
APP_DIR = ./programs
SRC_DIR = ./src

GEN_OBJ_DIR = $(OBJ_DIR)/general
GEN_APP_DIR = $(APP_DIR)/general
GEN_SRC_DIR = $(SRC_DIR)/general

GEN_SRC = $(wildcard $(GEN_SRC_DIR)/*.cpp)
GEN_OBJECTS = $(GEN_SRC:$(GEN_SRC_DIR)/%.cpp=$(GEN_OBJ_DIR)/%.o)


MKL_OBJ_DIR = $(OBJ_DIR)/mkl
MKL_APP_DIR =$(APP_DIR)/mkl
MKL_SRC_DIR = $(SRC_DIR)/mkl

MKL_SRC = $(wildcard $(MKL_SRC_DIR)/*.cpp)
MKL_OBJECTS = $(GEN_OBJECTS) $(MKL_SRC:$(MKL_SRC_DIR)/%.cpp=$(MKL_OBJ_DIR)/%.o)


CUDA_OBJ_DIR = $(OBJ_DIR)/cuda
CUDA_APP_DIR =$(APP_DIR)/cuda
CUDA_SRC_DIR = $(SRC_DIR)/cuda

CUDA_SRC = $(wildcard $(CUDA_SRC_DIR)/*.cpp)
CUDA_OBJECTS := $(GEN_OBJECTS) $(CUDA_SRC: $(CUDA_SRC_DIR)/%.cpp=$(CUDA_OBJ_DIR)/%.o)


all:  mkl_test cuda_test 


$(GEN_OBJ_DIR)/%.o : $(GEN_SRC_DIR)/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $@ -c $<	

$(MKL_OBJ_DIR)/%.o : $(MKL_SRC_DIR)/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(MKL_CXXFLAGS) $(MKL_INCLUDE) $(MKL_LIBRARY) -o $@ -c $<

$(MKL_APP_DIR)/$(MKL_TARGET) : $(MKL_OBJECTS)
	@mkdir -p $(@D)
	$(CXX) $(MKL_CXXFLAGS) $(MKL_INCLUDE) $(MKL_LDFLAGS) -o $(MKL_APP_DIR)/$(MKL_TARGET) $(MKL_OBJECTS)

$(CUDA_OBJ_DIR)/%.o : $(CUDA_SRC_DIR)/%.cpp
	@mkdir -p $(@D)
	$(NVCC) $(CUDA_CXXFLAGS) $(CUDA_INCLUDE) $(CUDA_LIBRARY) -o $@ -c $<

$(CUDA_APP_DIR)/$(CUDA_TARGET) : $(CUDA_OBJECTS)
	@mkdir -p $(@D)
	$(NVCC) $(CUDA_CXXFLAGS) $(CUDA_INCLUDE) $(CUDA_LDFLAGS) -o $(CUDA_APP_DIR)/$(CUDA_TARGET) $(CUDA_OBJECTS)





.PHONY: all build clean general mkl_build cuda_build

build: 
	@mkdir -p $(APP_DIR)
	@mkdir -p $(OBJ_DIR)
	@mkdir -p $(GEN_APP_DIR)
	@mkdir -p $(GEN_OBJ_DIR)

mkl_build: build $(MKL_APP_DIR)/$(MKL_TARGET)
	@mkdir -p $(MKL_APP_DIR)
	@mkdir -p $(MKL_OBJ_DIR)

cuda_build: build $(CUDA_APP_DIR)/$(CUDA_TARGET)
	@mkdir -p $(CUDA_APP_DIR)
	@mkdir -p $(CUDA_OBJ_DIR)

general: build $(GEN_OBJECTS)

mkl_test: $(MKL_APP_DIR)/$(MKL_TARGET)
	./$(MKL_APP_DIR)/$(MKL_TARGET)

cuda_test: $(CUDA_APP_DIR)/$(CUDA_TARGET)
	./$(CUDA_APP_DIR)/$(CUDA_TARGET)

clean:
	-@rm -rvf $(OBJ_DIR)/*
	-@rm -rvf $(APP_DIR)/*
