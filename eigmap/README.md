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

string                  reference eigenvector file
string                  reference frequency file
string                  eigenvector file to compare
string                  frequency file to compare
int                   number of atoms to remap
int int               atom_label in ref. structure -> atom_label in comp. structure
 ...

```

where `int` is an integer number and `string` a string of type char. The command line syntax can be shown by using the `-h` option:

```

$ ./eigmap -h
       _                               
   ___(_) __ _ _ __ ___   __ _ _ __    
  / _ \ |/ _` | '_ ` _ \ / _` | '_ \   
 |  __/ | (_| | | | | | | (_| | |_) |  
  \___|_|\__, |_| |_| |_|\__,_| .__/   
         |___/                |_|      
                              2.2

 Syntax: eigmap <setting file>

```

After the execution the files eigmap.txt and eigscal.txt are created. The header of the file contains information about the content.

## Example

The folder example contains the necessary files to test the code. The *eigmap.inp* file is an example of input file. The *.ascii* files may be used to visualize the phonon displacement pattern with, for example, the [v-sim](https://www.mem-lab.fr/en/Pages/L_SIM/Softwares/V_Sim.aspx) software.

## Citation

 The users of EIGMAP have little formal obligations specified in the [GNU General Public License](http://www.gnu.org/copyleft/gpl.txt).
 However, it is common practice in the scientific literature, to acknowledge the efforts of people that have made the research possible.
 In this spirit, please cite


A. Cammarata, M. Dasic and P. Nicolini, *Normal Dynamics: solving Newton’s equations in the reciprocal space*, Phys. Rev. Lett **XX**, XXXXX (XXXX) DOI: [xxx](https://doi.org/10.1103/xxx)

A. Cammarata, M. Dasic and P. Nicolini, *Sampling dynamical trajectories in the reciprocal space*, Phys. Rev. B **XX**, XXXXX (XXXX) DOI: [xxx](https://doi.org/10.1103/xxx)

