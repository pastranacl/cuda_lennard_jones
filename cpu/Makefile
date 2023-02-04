# Cesar L. Pastrana, 2021

CXX := g++
FLAGS := -std=c++11 -O3 -lm 
FLAGS_DEBUG := -g -Wall
OPENMP_FLAG := -fopenmp

BIN_PATH := ./bin
SRC_PATH := ./src

SRCS := $(wildcard ${SRC_PATH}/*.cpp ${SRC_PATH}/*.c)

ON := tubenz

parallel: 
	@clear
	@echo -n "Compiling parallel code... "
	@${CXX} ${SRCS} ${FLAGS} ${OPENMP_FLAG} -o ${BIN_PATH}/${ON}
	@echo "Done!"
	@echo "You can select the number of threads to use attending to the number of cores of the computer. For it, set the state variable export OMP_NUM_THREADS=<number of threads to use>"

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



