#!/bin/sh

plotShapelets.py -n 4 -p -b 1.0,1.0,0.44

solveShapelet.py ../data/N6251_test.fits -r 1028,1097,1025,1074 -N 967,1067,972,1026 -n 20 -m polar
solveShapelet.py ../data/N6251_test.fits -r 1028,1097,1025,1074 -N 967,1067,972,1026 -n 20 -b 3.,3. -p 0.448266


