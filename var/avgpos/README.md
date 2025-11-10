# avgpos

A tool to calculate the average position and standard deviation of selected atoms along a specified crystallographic direction from a POSCAR file.

## Features

- Read POSCAR files (VASP structure format)
- Select atoms by element type or by indices
- Calculate average position along:
  - Cartesian directions (x, y, z)
  - Crystallographic lattice vectors (a, b, c)
  - Custom Miller indices [h,k,l]
- Calculate standard deviation of positions
- Display individual atomic positions along the selected direction

## Requirements

- Python 3.6 or higher
- NumPy

## Installation

No installation required. Simply make the script executable:

```bash
chmod +x avgpos.py
```

Or run it with Python:

```bash
python3 avgpos.py
```

## Usage

### Basic syntax

```bash
./avgpos.py POSCAR -s <elements> -d <direction>
./avgpos.py POSCAR -i <indices> -d <direction>
```

### Options

- `POSCAR`: Path to the POSCAR file (required)
- `-s, --select`: Select atoms by element symbol(s), comma-separated (e.g., "Se" or "W,Mo")
- `-i, --indices`: Select atoms by indices (1-based), comma-separated (e.g., "1,2,3")
- `-d, --direction`: Direction specification (required):
  - Cartesian: `x`, `y`, `z`
  - Lattice vectors: `a`, `b`, `c`
  - Miller indices: `[h,k,l]` (e.g., `[1,1,0]`)

### Examples

Calculate average position of all Se atoms along the z-axis:
```bash
./avgpos.py POSCAR -s Se -d z
```

Calculate average position of atoms 2, 3, and 4 along the c lattice vector:
```bash
./avgpos.py POSCAR -i 2,3,4 -d c
```

Calculate average position of W and Mo atoms along the [1,1,0] direction:
```bash
./avgpos.py POSCAR -s W,Mo -d [1,1,0]
```

Calculate average position of all atoms of multiple elements along x-axis:
```bash
./avgpos.py POSCAR -s Se,Mo -d x
```

## Output

The tool displays:
- Structure information (number and types of atoms)
- Selected atoms
- Direction vector in Cartesian coordinates
- Average position along the direction (in Ångströms)
- Standard deviation (in Ångströms)
- Individual atomic positions (if 20 or fewer atoms are selected)

## Example Output

```
Reading POSCAR file: POSCAR
Structure contains 6 atoms:
  W: 1
  Se: 4
  Mo: 1

Selected 4 atom(s) of type: Se
Direction vector (Cartesian): [0.000000, 0.000000, 1.000000]

============================================================
RESULTS
============================================================
Number of atoms: 4
Average position: 35.123456 Å
Standard deviation: 2.345678 Å

Individual positions along direction:
  Atom 2: 38.343795 Å
  Atom 3: 28.534648 Å
  Atom 4: 31.870770 Å
  Atom 5: 35.012611 Å
```

## POSCAR File Format

The tool supports standard VASP POSCAR format with:
- Comment line
- Scale factor (preferably 1.0)
- Lattice vectors (3 lines)
- Element symbols
- Atom counts per element
- Coordinate type (Direct/Cartesian)
- Atomic positions

Both Direct (fractional) and Cartesian coordinates are supported.

## License

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

## Author

Part of the phtools collection: https://github.com/acammarat/phtools
