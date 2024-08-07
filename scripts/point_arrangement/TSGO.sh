#!/bin/bash

SRC=main.cpp
SRC1=CG_descent.cpp
SRC2=ObjGrad.cpp
SRC3=Variables.cpp
SRC4=Shake.cpp
SRC5=ConstructNeighbor.cpp
SRC6=LocalSearch.cpp
SRC7=MBH.cpp
SRC8=InitialSolution.cpp
SRC9=InOutput.cpp
SRC10=TSGO.cpp
SRC11=TabuSearch.cpp
EXE=TSGOP

srun g++ ${SRC} ${SRC1} ${SRC2} ${SRC3} ${SRC4} ${SRC5} ${SRC6} ${SRC7} ${SRC8} ${SRC9} ${SRC10} ${SRC11} -O3 -o ${EXE}


