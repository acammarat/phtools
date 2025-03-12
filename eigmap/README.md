# eigmap

Program to calculate the map between different eigenvector via scalar products

## Installation

The code requires a fortran compiler. After cloning, enter the folder and compile it with

`make`

If the compilation ends successfully, the executable eigmap is created.

## Usage

The qmatrix.nd and freq.nd files of the reference and of the structure to be compared must be be created with the tool [qpoints](https://github.com/acammarat/phtools/tree/main/qpoints)

The format of the input file is


```
 char                       reference POSCAR file
 int                        number of atomic types in reference
 char double                atom symbol, mass [amu]
 char                       reference eigenvector file
 char                       reference frequency file
 int int int int            acoustic modes in reference: q, j1, j2, j3
 char                       comparison POSCAR file
 int                        number of atomic types in comparison
 char double                atom symbol, mass [amu]
 char                       comparison eigenvector file
 char                       comparison frequency file
 int int int int            acoustic modes in comparison: q, j1, j2, j3
 int                        number of atoms to map
 int int                    atom_label in ref. structure -> atom_label in comp. structure
 ... ...


```

The command line syntax can be shown by using the `-h` option:

```

$ ./eigmap -h
       _                               
   ___(_) __ _ _ __ ___   __ _ _ __    
  / _ \ |/ _` | '_ ` _ \ / _` | '_ \   
 |  __/ | (_| | | | | | | (_| | |_) |  
  \___|_|\__, |_| |_| |_|\__,_| .__/   
         |___/                |_|      
                              2.3

 Syntax: eigmap <setting file>

```

After the execution the files eigmap.txt and eigscal.txt are created. The header of the file contains information about the content.

## Note

At the moment, the code works only for eigenvectors at Gamma.

## Example

The folder example contains the necessary files to test the code. The *eigmap.inp* file is an example of input file. The *.ascii* files may be used to visualize the phonon displacement pattern with, for example, the [v-sim](https://www.mem-lab.fr/en/Pages/L_SIM/Softwares/V_Sim.aspx) software.

## Citation

 The users of EIGMAP have little formal obligations specified in the [GNU General Public License](http://www.gnu.org/copyleft/gpl.txt).
 However, it is common practice in the scientific literature, to acknowledge the efforts of people that have made the research possible.
 In this spirit, please cite


A. Cammarata, M. Dasic and P. Nicolini, *Integrating Newton's equations of motion in the reciprocal space*, Journal of Chemical Physics, **161**, 084111 (2024) DOI: [10.1063/5.0224108](https://doi.org/10.1063/5.0224108)

where the formulation used to perform normal dynamics is reported.

