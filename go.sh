#!/bin/bash

#g++ main_1.cpp solver.cpp cglv.cpp assign_1.cpp wyarray.cpp inout.cpp
#g++ main_2.cpp readpara.cpp solverlin.cpp cglv.cpp cgproj.cpp pdmtrfunc.cpp kktecqfunc.cpp wyarray.cpp inout.cpp

sub="readpara.cpp wyarray.cpp inout.cpp"
g++ main_cg.cpp cglv.cpp pdmatrix.cpp $sub
