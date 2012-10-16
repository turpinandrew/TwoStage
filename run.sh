#!/bin/sh

#
# Run for various values of DECISION_RANK 4
#

for i in 1 2 3 4 5 6 7 8 9 10
do
    cp model.h model.h.temp
    sed 's/#define DECISION_RANK.*$/#define DECISION_RANK '"$i"'/' < model.h.temp > model.h
    make
    cat model.h
    ./model x
    ./model x
    ./model x
    ./model x
    ./model x
done
