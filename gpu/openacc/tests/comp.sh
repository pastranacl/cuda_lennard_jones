g++ main.cpp -O3 -fopenacc -fopenacc-dim=1024:1:256 -fcf-protection=none -foffload=-misa=sm_35 -o test
# -fopenacc-dim=1024:1:128 1024 gangs, 1 worker. 128 vectors
