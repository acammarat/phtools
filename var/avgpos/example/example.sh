#!/bin/bash
# Example usage of avgpos tool

echo "====== Example 1: Average position of Se atoms along z-axis ======"
python3 ../avgpos.py POSCAR -s Se -d z

echo ""
echo "====== Example 2: Average position of atoms 2,3,4 along c lattice vector ======"
python3 ../avgpos.py POSCAR -i 2,3,4 -d c

echo ""
echo "====== Example 3: Average position of W and Mo atoms along [1,1,0] direction ======"
python3 ../avgpos.py POSCAR -s W,Mo -d "[1,1,0]"

echo ""
echo "====== Example 4: Average position of all Se atoms along x-axis ======"
python3 ../avgpos.py POSCAR -s Se -d x
