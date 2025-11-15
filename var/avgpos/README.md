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
- **NEW:** Calculate and export orthogonal projections onto a plane perpendicular to the direction vector

## Requirements

- Python 3.6 or higher
- NumPy
- Matplotlib (for generating heatmap visualizations)
- SciPy (for interpolation in smooth heatmaps)

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
- `-o, --output`: Output file for plane projection data (optional)
- `--plot`: Generate Python matplotlib script for heatmap visualization (requires `-o`)
- `--labels`: Include atom labels (element+ID, e.g., Se2, Ti4) in output and plot (requires `-o`)
- `--replicate`: Replicate the plot along e and f axes (format: "ne,nf", default: "1,1")
  - Supports non-integer replication (e.g., "2.5,3" for 2.5x3 replication)
- `--no-circles`: Hide circles representing atom positions in the plot (only effective when `--labels` is not used)

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

Calculate average position and export plane projection data:
```bash
./avgpos.py POSCAR -s Se -d z -o projections.dat
```

Calculate average position and generate matplotlib heatmap script:
```bash
./avgpos.py POSCAR -s Se -d z -o projections.dat --plot
# Then run: python3 projections_plot.py
```

Calculate average position with atom labels and generate labeled heatmap:
```bash
./avgpos.py POSCAR -s Se -d z -o projections.dat --plot --labels
# Then run: python3 projections_plot.py
# Labels will show atom type and POSCAR file ID (e.g., Se2, Se3, Ti1)
```

Generate heatmap with 2.5x3 replication along e and f axes:
```bash
./avgpos.py POSCAR -s Se -d z -o projections.dat --plot --replicate 2.5,3
# Then run: python3 projections_plot.py
# Plot will show 2.5 replications along e-axis and 3 along f-axis
```

Generate smooth heatmap without atom position circles:
```bash
./avgpos.py POSCAR -s Se -d z -o projections.dat --plot --no-circles
# Then run: python3 projections_plot.py
# Plot will show only the smooth interpolated surface without circles
```

## Output

### Standard Output

The tool displays:
- Structure information (number and types of atoms)
- Selected atoms
- Direction vector in Cartesian coordinates
- Average position along the direction (in Ångströms)
- Standard deviation (in Ångströms)
- Individual atomic positions (if 20 or fewer atoms are selected)

### Plane Projection Output File (optional)

When the `-o` option is specified, the tool generates a data file containing:
- **Column 1 (e)**: First coordinate of the atom's orthogonal projection onto the plane
- **Column 2 (f)**: Second coordinate of the atom's orthogonal projection onto the plane  
- **Column 3 (g)**: Signed distance from the plane (average_position - atom_distance_along_direction)
- **Column 4 (label)**: Atom label with element type and unique ID (e.g., Ti1, O2) - only when `--labels` is used

The plane is perpendicular to the specified direction vector and passes through the calculated average position. The e and f coordinates form an orthonormal 2D coordinate system in the plane.

### Matplotlib Plotting Script (optional)

When the `--plot` flag is used along with `-o`, the tool generates a Python script using matplotlib that creates a smooth interpolated heatmap visualization of the plane projection data:
- **Script file**: Named as `<output_basename>_plot.py`
- **Output image**: Named as `<output_basename>_heatmap.png`
- The heatmap uses Radial Basis Function (RBF) interpolation to create a smooth surface covering the entire e,f range
- g values are represented by a color gradient (coolwarm colormap)
- Original data points are overlaid as black dots for reference
- When `--labels` is also used, atom labels (element+POSCAR file ID) are annotated on the plot
- To generate the plot, run: `python3 <script_file>`

**Requirements**: Matplotlib and SciPy must be installed on your system to generate the visualization.

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
