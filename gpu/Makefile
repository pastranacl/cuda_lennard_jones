# Cesar L. Pastrana, 2021


CXX := nvcc


BIN_PATH := ./bin
SRC_PATH := ./src

SRCS := $(wildcard ${SRC_PATH}/*.cpp wildcard ${SRC_PATH}/*.cu)

ON := gputubenz

main: 
	@echo "${SRCS}"
	@clear
	@echo -n "Compiling GPU code... "
	@${CXX} ${SRCS} -o ${BIN_PATH}/${ON}
	@echo -n "Done!\n"
	



