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


def get_atom_labels(structure, atom_indices):
    """
    Generate atom labels with element type and POSCAR file ID for selected atoms.
    
    Parameters:
    -----------
    structure : dict
        Structure dictionary from read_poscar
    atom_indices : numpy.ndarray
        Indices of atoms (0-based)
        
    Returns:
    --------
    list : List of atom labels (e.g., ['Ti1', 'Ti2', 'O3', 'O4'])
        where the number corresponds to the atom's position in the POSCAR file (1-based)
    """
    labels = []
    
    # Create a mapping of atom index to element type
    idx = 0
    atom_to_element = []
    for i, element in enumerate(structure['elements']):
        count = structure['atom_counts'][i]
        atom_to_element.extend([element] * count)
        idx += count
    
    # Use the POSCAR file index (1-based) as the atom ID
    for atom_idx in atom_indices:
        element = atom_to_element[atom_idx]
        # atom_idx is 0-based, so add 1 to get POSCAR file position
        labels.append(f"{element}{atom_idx + 1}")
    
    return labels


def generate_plot_script(data_file, script_file, output_image='heatmap.png', with_labels=False, replicate=(1, 1), no_circles=False):
    """
    Generate a Python script using matplotlib to plot the plane projection data as a heatmap.
    
    Parameters:
    -----------
    data_file : str
        Path to the data file containing projection data
    script_file : str
        Path where the Python plotting script will be written
    output_image : str
        Name of the output image file (default: 'heatmap.png')
    with_labels : bool
        Whether to include atom labels in the plot
    replicate : tuple
        Number of replications along e and f axes (ne, nf)
    no_circles : bool
        Whether to hide the circles representing atom positions (only when labels are not used)
    """
    script_content = f"""#!/usr/bin/env python3
\"\"\"
Matplotlib script to visualize plane projection data as a smooth interpolated heatmap.
Generated automatically by avgpos tool.

Usage: python3 {script_file}
\"\"\"

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.interpolate import Rbf
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Read data from file
data = np.loadtxt('{data_file}', dtype=str)

# Extract coordinates and g values
e_orig = data[:, 0].astype(float)
f_orig = data[:, 1].astype(float)
g_orig = data[:, 2].astype(float)

# Extract labels if present (4th column)
has_labels = data.shape[1] > 3
if has_labels:
    labels_orig = data[:, 3]

# Replication parameters
ne_rep, nf_rep = {replicate[0]}, {replicate[1]}

# Replicate data along e and f axes
e_list, f_list, g_list, labels_list = [], [], [], []
e_range_orig = e_orig.max() - e_orig.min() if len(e_orig) > 1 else 1.0
f_range_orig = f_orig.max() - f_orig.min() if len(f_orig) > 1 else 1.0

# Determine how many full and partial replications to make
ne_full = int(np.ceil(ne_rep))
nf_full = int(np.ceil(nf_rep))

for ie in range(ne_full):
    for jf in range(nf_full):
        # Calculate if this replica is fully or partially included
        e_factor = min(1.0, ne_rep - ie) if ie < ne_full - 1 else (ne_rep - ie)
        f_factor = min(1.0, nf_rep - jf) if jf < nf_full - 1 else (nf_rep - jf)
        
        # Include this replica if it has non-zero contribution
        if e_factor > 0 and f_factor > 0:
            e_shift = ie * (e_range_orig + (e_orig.max() - e_orig.min()))
            f_shift = jf * (f_range_orig + (f_orig.max() - f_orig.min()))
            e_list.append(e_orig + e_shift)
            f_list.append(f_orig + f_shift)
            g_list.append(g_orig)
            if has_labels:
                labels_list.append(labels_orig)

e = np.concatenate(e_list)
f = np.concatenate(f_list)
g = np.concatenate(g_list)
if has_labels:
    labels = np.concatenate(labels_list)

# Create a regular grid for interpolation
# Determine the range of e and f
e_min, e_max = e.min(), e.max()
f_min, f_max = f.min(), f.max()

e_grid = np.linspace(e_min, e_max, 200)
f_grid = np.linspace(f_min, f_max, 200)
e_mesh, f_mesh = np.meshgrid(e_grid, f_grid)

# Use Radial Basis Function interpolation with very small smoothing
# This ensures interpolation passes extremely close to data points while handling duplicates
# smooth value is set to a tiny value to get nearly exact values at data points
try:
    rbf = Rbf(e, f, g, function='thin_plate', smooth=1e-10)
    g_interp = rbf(e_mesh, f_mesh)
except:
    # Fall back to multiquadric if thin_plate fails
    try:
        rbf = Rbf(e, f, g, function='multiquadric', smooth=1e-10)
        g_interp = rbf(e_mesh, f_mesh)
    except:
        # Last resort: use linear with small smoothing
        rbf = Rbf(e, f, g, function='linear', smooth=1e-8)
        g_interp = rbf(e_mesh, f_mesh)

# Determine color range from actual data values (not interpolated)
vmin, vmax = g.min(), g.max()

# Create figure and axis
fig, ax = plt.subplots(figsize=(10, 8))

# Create smooth heatmap using pcolormesh with RGB gradient (jet colormap)
# Use the same vmin/vmax as the scatter plot for consistent colors
heatmap = ax.pcolormesh(e_mesh, f_mesh, g_interp, cmap='jet', shading='auto', vmin=vmin, vmax=vmax)
"""
    
    # Add scatter points conditionally
    if with_labels or not no_circles:
        script_content += f"""
# Overlay the original data points with their EXACT g values colored
# This ensures atomic positions correspond to the real g value from the data file
scatter = ax.scatter(e, f, c=g, cmap='jet', s=150, edgecolors='black', linewidths=2, zorder=10, vmin=vmin, vmax=vmax)
"""
    
    script_content += f"""
# Add colorbar with height matching the plot
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.1)
cbar = plt.colorbar(heatmap, cax=cax)
cbar.set_label('Distance from plane (Å)', fontsize=12)

# Set axis labels with units
ax.set_xlabel('x (Å)', fontsize=12)
ax.set_ylabel('y (Å)', fontsize=12)
"""
    
    if with_labels:
        script_content += f"""
# Add atom labels (already replicated above if needed)
if has_labels:
    for i in range(len(e)):
        ax.annotate(labels[i], (e[i], f[i]), 
                    xytext=(5, 5), textcoords='offset points',
                    fontsize=10, fontweight='bold', color='black',
                    bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='black', alpha=0.7))
"""
    
    script_content += f"""
# Add grid
ax.grid(True, alpha=0.3)

# Set equal aspect ratio
ax.set_aspect('equal', adjustable='box')

# Save figure
plt.tight_layout()
plt.savefig('{output_image}', dpi=150, bbox_inches='tight')
print(f"Plot saved to {output_image}")

# Optionally display the plot (comment out if running headless)
# plt.show()
"""
    
    with open(script_file, 'w') as f:
        f.write(script_content)
    
    # Make the script executable (Unix-like systems)
    import os
    try:
        os.chmod(script_file, 0o755)
    except:
        pass  # Ignore if chmod fails (e.g., on Windows)


def calculate_plane_projections(structure, atom_indices, direction_vector, average_position):
    """
    Calculate orthogonal projections of atoms onto a plane perpendicular to 
    the direction vector and passing through the average position.
    
    Parameters:
    -----------
    structure : dict
        Structure dictionary from read_poscar
    atom_indices : numpy.ndarray
        Indices of atoms to include in calculation
    direction_vector : numpy.ndarray
        Unit vector defining the direction (normal to the plane)
    average_position : float
        Average position along the direction (defines plane location)
        
    Returns:
    --------
    numpy.ndarray : Nx3 array where each row contains [e, f, g]
        - e, f: 2D coordinates of the projection on the plane
        - g: average_position minus the distance of the atom from the plane
    """
    positions = structure['positions'][atom_indices]
    
    # Calculate the distance of each atom along the direction vector
    distances_along_dir = np.dot(positions, direction_vector)
    
    # Calculate the signed distance from each atom to the plane
    # (positive if atom is on the side of the direction vector, negative otherwise)
    signed_distances = distances_along_dir - average_position
    
    # Project each atom onto the plane
    # projection = position - (signed_distance * normal_vector)
    projections_3d = positions - np.outer(signed_distances, direction_vector)
    
    # Create an orthonormal basis for the plane
    # Find two orthogonal vectors in the plane
    # Start with an arbitrary vector not parallel to direction_vector
    if abs(direction_vector[2]) < 0.9:
        arbitrary = np.array([0.0, 0.0, 1.0])
    else:
        arbitrary = np.array([1.0, 0.0, 0.0])
    
    # First basis vector in the plane (orthogonal to direction_vector)
    basis1 = arbitrary - np.dot(arbitrary, direction_vector) * direction_vector
    basis1 = basis1 / np.linalg.norm(basis1)
    
    # Second basis vector (orthogonal to both direction_vector and basis1)
    basis2 = np.cross(direction_vector, basis1)
    basis2 = basis2 / np.linalg.norm(basis2)
    
    # Project the 3D projections onto the 2D plane coordinate system
    e_coords = np.dot(projections_3d, basis1)
    f_coords = np.dot(projections_3d, basis2)
    
    # Calculate g = average_position - distance_from_plane
    # The distance from the plane is the absolute value of signed_distances
    # But we want: average_position - distance_of_atom_from_plane
    # Since distance_along_dir = average_position + signed_distance
    # we have: g = average_position - abs(signed_distance) if we want actual distance
    # But the requirement says: "average_position minus the distance of the atom from the plane"
    # which could mean: average_position - distance_along_dir = -signed_distances
    g_coords = -signed_distances
    
    # Combine into Nx3 array
    result = np.column_stack((e_coords, f_coords, g_coords))
    
    return result


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
    parser.add_argument('-o', '--output', type=str,
                        help='Output file for plane projection data (3 columns: e, f, g)')
    parser.add_argument('--plot', action='store_true',
                        help='Generate Python matplotlib script for heatmap visualization (requires -o)')
    parser.add_argument('--labels', action='store_true',
                        help='Include atom labels (element+ID) in output and plot (requires -o)')
    parser.add_argument('--replicate', type=str, default='1,1',
                        help='Replicate the plot along e and f axes (format: "ne,nf", e.g., "2.5,3" for 2.5x3 replication)')
    parser.add_argument('--no-circles', action='store_true',
                        help='Hide circles representing atom positions in the plot (only when --labels is not used)')
    
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
    
    # Calculate and write plane projections if output file is specified
    if args.output:
        projections = calculate_plane_projections(
            structure, atom_indices, direction_vector, average
        )
        
        # Get atom labels if requested
        if args.labels:
            labels = get_atom_labels(structure, atom_indices)
            # Create a combined array with projections and labels
            # Save with labels as 4th column
            with open(args.output, 'w') as f:
                f.write('# e f g label\n')
                f.write('# Projections onto plane perpendicular to direction vector\n')
                f.write('# e, f: 2D coordinates on plane\n')
                f.write('# g: average_position - distance_from_plane\n')
                f.write('# label: atom type and ID (e.g., Ti1, O2)\n')
                for i, label in enumerate(labels):
                    f.write(f"{projections[i, 0]:.6f} {projections[i, 1]:.6f} {projections[i, 2]:.6f} {label}\n")
            
            print()
            print(f"Plane projection data with labels written to: {args.output}")
            print(f"  Columns: e, f, g, label")
        else:
            # Write to file without labels
            np.savetxt(args.output, projections, fmt='%.6f', 
                       header='e f g\nProjections onto plane perpendicular to direction vector\n'
                              'e, f: 2D coordinates on plane\n'
                              'g: average_position - distance_from_plane',
                       comments='# ')
            
            print()
            print(f"Plane projection data written to: {args.output}")
            print(f"  Columns: e, f, g")
        
        print(f"  e, f: 2D coordinates of atom projection on plane")
        print(f"  g: signed distance from plane (average_position - atom_distance)")
        
        # Generate matplotlib plot script if requested
        if args.plot:
            import os
            # Determine output paths
            base_name = os.path.splitext(args.output)[0]
            script_file = f"{base_name}_plot.py"
            image_file = f"{base_name}_heatmap.png"
            
            # Parse replication argument
            try:
                replicate_parts = args.replicate.split(',')
                ne_rep = float(replicate_parts[0])
                nf_rep = float(replicate_parts[1]) if len(replicate_parts) > 1 else ne_rep
                replicate = (ne_rep, nf_rep)
            except (ValueError, IndexError):
                print(f"Warning: Invalid replicate format '{args.replicate}'. Using default 1,1")
                replicate = (1, 1)
            
            # Generate the plotting script
            generate_plot_script(args.output, script_file, image_file, args.labels, replicate, args.no_circles)
            
            print()
            print(f"Matplotlib plotting script generated: {script_file}")
            print(f"To create the heatmap, run: python3 {script_file}")
            print(f"Output image will be: {image_file}")
            if args.labels:
                print(f"  (with atom labels)")
            if replicate != (1, 1):
                print(f"  (with {replicate[0]}x{replicate[1]} replication)")
            if args.no_circles and not args.labels:
                print(f"  (without atom position circles)")
    elif args.plot or args.labels:
        print()
        print("Warning: --plot and --labels flags require -o/--output to be specified. Ignoring.")
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
