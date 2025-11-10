#!/usr/bin/env python3
"""
avgpos - Calculate average position and standard deviation of selected atoms
along a crystallographic direction from a POSCAR file.

Copyright (C) 2025
This file is part of phtools.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
"""

import sys
import numpy as np
import argparse


def read_poscar(filename):
    """
    Read a POSCAR file and return atomic structure information.
    
    Parameters:
    -----------
    filename : str
        Path to the POSCAR file
        
    Returns:
    --------
    dict : Dictionary containing structure information
        - 'lattice': 3x3 numpy array with lattice vectors
        - 'elements': list of element symbols
        - 'atom_counts': list of atom counts per element
        - 'positions': Nx3 numpy array with atomic positions
        - 'coordinate_type': 'Direct' or 'Cartesian'
    """
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    # Read comment line
    comment = lines[0].strip()
    
    # Read scale factor
    scale = float(lines[1].strip())
    if abs(scale - 1.0) > 1e-6:
        print(f"Warning: Scale factor is {scale}, not 1.0. Applying scaling.")
    
    # Read lattice vectors
    lattice = np.zeros((3, 3))
    for i in range(3):
        lattice[i] = [float(x) for x in lines[2+i].split()]
    lattice *= scale
    
    # Read element symbols
    elements = lines[5].split()
    
    # Read atom counts
    atom_counts = [int(x) for x in lines[6].split()]
    total_atoms = sum(atom_counts)
    
    # Check for selective dynamics
    line_idx = 7
    if lines[line_idx].strip()[0].upper() in ['S']:
        line_idx += 1
    
    # Read coordinate type
    coord_type = lines[line_idx].strip()
    coordinate_type = 'Direct' if coord_type[0].upper() in ['D'] else 'Cartesian'
    
    # Read atomic positions
    positions = np.zeros((total_atoms, 3))
    for i in range(total_atoms):
        pos_line = lines[line_idx + 1 + i].split()
        positions[i] = [float(x) for x in pos_line[:3]]
    
    # Convert direct coordinates to Cartesian if needed
    if coordinate_type == 'Direct':
        positions = np.dot(positions, lattice)
    
    return {
        'lattice': lattice,
        'elements': elements,
        'atom_counts': atom_counts,
        'positions': positions,
        'coordinate_type': coordinate_type
    }


def select_atoms(structure, selection):
    """
    Select atoms based on element type or indices.
    
    Parameters:
    -----------
    structure : dict
        Structure dictionary from read_poscar
    selection : str or list
        Either element symbol(s) or atom indices (1-based)
        
    Returns:
    --------
    numpy.ndarray : Indices of selected atoms (0-based)
    """
    if isinstance(selection, str):
        # Selection by element
        selected_indices = []
        idx = 0
        for i, element in enumerate(structure['elements']):
            count = structure['atom_counts'][i]
            if element in selection.split(','):
                selected_indices.extend(range(idx, idx + count))
            idx += count
        return np.array(selected_indices)
    else:
        # Selection by indices (convert from 1-based to 0-based)
        return np.array([i-1 for i in selection])


def get_direction_vector(structure, direction):
    """
    Get the unit vector for the specified crystallographic direction.
    
    Parameters:
    -----------
    structure : dict
        Structure dictionary from read_poscar
    direction : str
        Direction specification: 'x', 'y', 'z', 'a', 'b', 'c', or custom [h,k,l]
        
    Returns:
    --------
    numpy.ndarray : Unit vector in Cartesian coordinates
    """
    lattice = structure['lattice']
    
    if direction.lower() == 'x':
        return np.array([1.0, 0.0, 0.0])
    elif direction.lower() == 'y':
        return np.array([0.0, 1.0, 0.0])
    elif direction.lower() == 'z':
        return np.array([0.0, 0.0, 1.0])
    elif direction.lower() == 'a':
        vec = lattice[0]
    elif direction.lower() == 'b':
        vec = lattice[1]
    elif direction.lower() == 'c':
        vec = lattice[2]
    else:
        # Parse custom direction [h,k,l]
        try:
            direction = direction.strip('[]')
            h, k, l = [float(x) for x in direction.split(',')]
            # Convert Miller indices to Cartesian
            vec = h * lattice[0] + k * lattice[1] + l * lattice[2]
        except:
            raise ValueError(f"Invalid direction specification: {direction}")
    
    # Normalize to unit vector
    return vec / np.linalg.norm(vec)


def calculate_average_position(structure, atom_indices, direction_vector):
    """
    Calculate average position and standard deviation along a direction.
    
    Parameters:
    -----------
    structure : dict
        Structure dictionary from read_poscar
    atom_indices : numpy.ndarray
        Indices of atoms to include in calculation
    direction_vector : numpy.ndarray
        Unit vector defining the direction
        
    Returns:
    --------
    tuple : (average, std_dev, positions_along_dir)
        - average: mean position along the direction
        - std_dev: standard deviation
        - positions_along_dir: array of positions projected onto direction
    """
    positions = structure['positions'][atom_indices]
    
    # Project positions onto the direction vector
    positions_along_dir = np.dot(positions, direction_vector)
    
    # Calculate statistics
    average = np.mean(positions_along_dir)
    std_dev = np.std(positions_along_dir)
    
    return average, std_dev, positions_along_dir


def main():
    parser = argparse.ArgumentParser(
        description='Calculate average position and standard deviation of selected atoms '
                    'along a crystallographic direction from a POSCAR file.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Average position of all Se atoms along z-axis
  %(prog)s POSCAR -s Se -d z
  
  # Average position of atoms 2,3,4 along c lattice vector
  %(prog)s POSCAR -i 2,3,4 -d c
  
  # Average position along custom direction [1,1,0]
  %(prog)s POSCAR -s W,Mo -d [1,1,0]
        """
    )
    
    parser.add_argument('poscar', help='Path to POSCAR file')
    parser.add_argument('-s', '--select', type=str, 
                        help='Select atoms by element symbol(s), comma-separated (e.g., "Se" or "W,Mo")')
    parser.add_argument('-i', '--indices', type=str,
                        help='Select atoms by indices (1-based), comma-separated (e.g., "1,2,3")')
    parser.add_argument('-d', '--direction', type=str, required=True,
                        help='Direction: x, y, z (Cartesian) or a, b, c (lattice vectors) '
                             'or [h,k,l] (Miller indices)')
    
    args = parser.parse_args()
    
    # Validate input
    if not args.select and not args.indices:
        parser.error("Must specify either --select or --indices")
    if args.select and args.indices:
        parser.error("Cannot specify both --select and --indices")
    
    # Read POSCAR file
    print(f"Reading POSCAR file: {args.poscar}")
    try:
        structure = read_poscar(args.poscar)
    except Exception as e:
        print(f"Error reading POSCAR file: {e}")
        sys.exit(1)
    
    # Print structure information
    total_atoms = sum(structure['atom_counts'])
    print(f"Structure contains {total_atoms} atoms:")
    for elem, count in zip(structure['elements'], structure['atom_counts']):
        print(f"  {elem}: {count}")
    print()
    
    # Select atoms
    if args.select:
        atom_indices = select_atoms(structure, args.select)
        print(f"Selected {len(atom_indices)} atom(s) of type: {args.select}")
    else:
        indices = [int(x) for x in args.indices.split(',')]
        atom_indices = select_atoms(structure, indices)
        print(f"Selected {len(atom_indices)} atom(s) by indices: {args.indices}")
    
    if len(atom_indices) == 0:
        print("Error: No atoms selected!")
        sys.exit(1)
    
    # Get direction vector
    try:
        direction_vector = get_direction_vector(structure, args.direction)
        print(f"Direction vector (Cartesian): [{direction_vector[0]:.6f}, "
              f"{direction_vector[1]:.6f}, {direction_vector[2]:.6f}]")
    except Exception as e:
        print(f"Error parsing direction: {e}")
        sys.exit(1)
    
    # Calculate average position
    average, std_dev, positions = calculate_average_position(
        structure, atom_indices, direction_vector
    )
    
    # Print results
    print()
    print("=" * 60)
    print("RESULTS")
    print("=" * 60)
    print(f"Number of atoms: {len(atom_indices)}")
    print(f"Average position: {average:.6f} Å")
    print(f"Standard deviation: {std_dev:.6f} Å")
    print()
    
    # Print individual positions if not too many
    if len(atom_indices) <= 20:
        print("Individual positions along direction:")
        for i, (idx, pos) in enumerate(zip(atom_indices, positions)):
            print(f"  Atom {idx+1}: {pos:.6f} Å")
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
