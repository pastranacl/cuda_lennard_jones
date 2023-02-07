# Cesar L. Pastrana, 2021

CXX := g++
FLAGS := -std=c++11 -O3 -lm 
FLAGS_DEBUG := -g -Wall
OPENMP_FLAGS := -fopenmp
OPEN_ACC_FLAGS := -fopenacc -fcf-protection=none -foffload=-misa=sm_35
#-fopenacc-dim=1024:1:256

BIN_PATH := ./bin/
SRC_PATH := ./src/

SRCS := $(wildcard ${SRC_PATH}/*.cpp ${SRC_PATH}/*.c)

ON := gputubenz

parallel: 
	@clear
	@echo -n "Compiling parallel code... "
	@${CXX} ${SRCS} ${FLAGS} ${OPEN_ACC_FLAGS} -o ${BIN_PATH}/${ON}
	@echo "Done!"

sequential:
	@clear
	@echo -n "Compiling sequential code... "
	@${CXX} ${SRCS} ${FLAGS} -o ${BIN_PATH}/${ON}
	@echo "Done!"

debug:
	@clear
	@echo -n "Compiling sequential code IN DEBUG MODE... "
	@${CXX} ${SRCS} ${FLAGS_DEBUG} -o ${BIN_PATH}/${ON}
	@echo "Done!"


install: filesystem parallel

filesystem:
# Create necessary folders 
ifeq ("$(wildcard ${BIN_PATH})","")
	@mkdir ${BIN_PATH}
endif

ifeq ("$(wildcard ${BIN_PATH}/ENZYMES)","")
	@mkdir ${BIN_PATH}/ENZYMES
endif


clean:
# Remove previous data
ifneq ("$(wildcard ${BIN_PATH}/ENZYMES/*.dat)","")
	@rm ${BIN_PATH}/ENZYMES/*.xyz
endif	



