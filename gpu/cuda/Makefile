# Cesar L. Pastrana, 2021


CXX := nvcc

CUFLAGS := -O3
BIN_PATH := ./bin
SRC_PATH := ./src

SRCS := $(wildcard ${SRC_PATH}/*.cpp wildcard ${SRC_PATH}/*.cu)

ON := gputubenz

main: 
	@echo "${SRCS}"
	@clear
	@echo -n "Compiling GPU code... "
	@${CXX} ${SRCS} ${CUFLAGS} -arch=sm_61 -o ${BIN_PATH}/${ON}
	@echo -n "Done!\n"
	



