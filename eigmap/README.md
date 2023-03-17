# eigmap

Program to calculate the map between different eigenvector via scalar products

## Installation

The code requires a fortran compiler. After cloning, enter the folder and compile it with

`make`

If the compilation ends successfully, the executable eigmap is created.

## Usage

The qmatrix.nd and freq.nd files of the reference and of the structure to be compared must be be created with the tool qpoints

`EIGENVECTORS = .TRUE.`

in the PHONOPY configuration file.

The syntax of the input file is


```

 int                     number of atomic groups
 int                     number of atoms in group 1
 int                     label of the first atom in the group
 ...
 int                     label of the last atom in the group
 int                     number of atoms in group 2
 int                     label of the first atom in the group
 ...
 int                     label of the last atom in the group
 ...                     number of atoms in group 3 - if any, and similar
                         input structure as above

```

where `int` is an integer number. The syntax can be shown by using the `-h` option:

```

$ ./phonchar -h
         _                      _
   _ __ | |__   ___  _ __   ___| |__   __ _ _ __
  | '_ \| '_ \ / _ \| '_ \ / __| '_ \ / _` | '__|
  | |_) | | | | (_) | | | | (__| | | | (_| | |
  | .__/|_| |_|\___/|_| |_|\___|_| |_|\__,_|_|
  |_|
                  2.2

 Syntax: phonchar <setting file> <band.yaml file>

```

After the execution the file phchar.dat is created. The header of the file contains information about the content.

## Example

The folder MoWS4 contains the necessary files to test the code. The *band.yaml* has been generated for the MoS2/WS2 bilayer structure reported in the *POSCAR* file; the *phchar.inp* file is an example of input file.

You can visualize the results by using [gnuplot](http://www.gnuplot.info/):

```

gnuplot> plot "phchar.dat" u 1:2:3 w p palette ps 0.5

```

The weight is rendered with the thermal map; the meaning of the weight value is explained in the header of the *phchar.dat* file.

If you run the example, you should get a figure like this

![example.png](example.png)

## Citation

 The users of PHONCHAR have little formal obligations specified in the [GNU General Public License](http://www.gnu.org/copyleft/gpl.txt).
 However, it is common practice in the scientific literature, to acknowledge the efforts of people that have made the research possible.
 In this spirit, please cite

 A. Cammarata, T. Polcar, *Fine control of lattice thermal conductivity in low-dimensional materials*, Phys. Rev. B **103**, 035406 (2021), DOI: [10.1103/PhysRevB.103.035406](https://doi.org/10.1103/PhysRevB.103.035406)

 where the formulation used to calculate the phonon atomic character is reported in section V "Atomic character of the phonon modes" of the Supplemental Material.

