#!/bin/bash
#LD_LIBRARY_PATH=clang/lib time ./sched ./instances/PSMF/E1/M3_N6_U1* | tee output.log
screen -L -S sched ./sched $*

