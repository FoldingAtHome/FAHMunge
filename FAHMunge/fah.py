##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2013 Stanford University and the Authors
#
# Authors: Kyle A. Beauchamp
# Contributors:
#
# MDTraj is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation, either version 2.1
# of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with MDTraj. If not, see <http://www.gnu.org/licenses/>.
##############################################################################

"""
Code for merging and munging trajectories from FAH datasets.
"""
##############################################################################
# imports
##############################################################################

from __future__ import print_function, division
import os
import glob
import tarfile
from mdtraj.formats.hdf5 import HDF5TrajectoryFile
import mdtraj as md
import tables
from mdtraj.utils.contextmanagers import enter_temp_directory
from mdtraj.utils import six

def keynat(string):
    '''A natural sort helper function for sort() and sorted()
    without using regular expression.

    >>> items = ('Z', 'a', '10', '1', '9')
    >>> sorted(items)
    ['1', '10', '9', 'Z', 'a']
    >>> sorted(items, key=keynat)
    ['1', '9', '10', 'Z', 'a']
    '''
    r = []
    for c in string:
        try:
            c = int(c)
            try:
                r[-1] = r[-1] * 10 + c
            except:
                r.append(c)
        except:
            r.append(c)
    return r

##############################################################################
# globals
##############################################################################


def strip_water(allatom_filename, protein_filename, protein_atom_indices, min_num_frames=1):
    """Strip water (or other) atoms from a Core17, Core18, or OCore FAH HDF5 trajectory.
    
    Parameters
    ----------
    allatom_filename : str
        Path to HDF5 trajectory with all atoms.  This trajectory must have been generated by
        concatenate_core17 or concatenate_siegetank--e.g. it must include
        extra metadata that lists the XTC files (bzipped or in OCore directories) that
        have already been processed.  This file will not be modified.
    protein_filename : str
        Path to HDF5 trajectory with all just protein atoms.  This trajectory must have been generated by
        concatenate_core17 or concatenate_siegetank--e.g. it must include
        extra metadata that lists the XTC files (bzipped or in OCore directories) that
        have already been processed.  This file will be appended to.
    protein_atom_indices : np.ndarray, dtype='int'
        List of atom indices to extract from allatom HDF5 file.
    min_num_frames : int, optional, default=1
        Skip if below this number.

    """    
    if not os.path.exists(allatom_filename):
        print("Skipping, %s not found" % allatom_filename)
        return

    trj_allatom = HDF5TrajectoryFile(allatom_filename, mode='r')
    
    print('all-atom trajectory %s has %d frames' % (allatom_filename, len(trj_allatom))) 
    if len(trj_allatom) < min_num_frames:
        print("Must have at least %d frames in %s to proceed!" % (min_num_frames, allatom_filename))
        del trj_allatom
        return

    if hasattr(trj_allatom.root, "processed_filenames"):
        key = "processed_filenames"  # Core17, Core18 style data
    elif hasattr(trj_allatom.root, "processed_directories"):
        key = "processed_directories"  # Siegetank style data
    else:
        raise(ValueError("Can't find processed files in %s" % allatom_filename))

    trj_protein = HDF5TrajectoryFile(protein_filename, mode='a')

    try:
        trj_protein._create_earray(where='/', name=key, atom=tables.StringAtom(1024), shape=(0,))
        trj_protein.topology = trj_allatom.topology.subset(protein_atom_indices)
    except tables.NodeError:
        pass

    n_frames_allatom = len(trj_allatom)
    try:
        n_frames_protein = len(trj_protein)
    except tables.NoSuchNodeError:
        n_frames_protein = 0

    filenames_allatom = getattr(trj_allatom.root, key)
    filenames_protein = getattr(trj_protein._handle.root, key)  # Hacky workaround of MDTraj bug #588
    
    n_files_allatom = len(filenames_allatom)
    n_files_protein = len(filenames_protein)
    print("Found %d,%d filenames and %d,%d frames in %s and %s, respectively." % (n_files_allatom, n_files_protein, n_frames_allatom, n_frames_protein, allatom_filename, protein_filename))
    
    if n_frames_protein > n_frames_allatom:
        raise(ValueError("Found more frames in protein trajectory (%d) than allatom trajectory (%d)" % (n_frames_protein, n_frames_allatom)))
    
    if n_files_protein > n_files_allatom:
        raise(ValueError("Found more filenames in protein trajectory (%d) than allatom trajectory (%d)" % (n_files_protein, n_files_allatom)))
    
    if n_frames_protein == n_frames_allatom or n_files_allatom == n_files_protein:
        if not (n_frames_protein == n_frames_allatom and n_files_allatom == n_files_protein):
            raise(ValueError("The trajectories must match in BOTH n_frames and n_filenames or NEITHER."))
        else:
            print("Same number of frames and filenames found, skipping.")
            del trj_allatom, trj_protein
            return

    trj_allatom.seek(n_frames_protein)  # Jump forward past what we've already stripped.
    coordinates, time, cell_lengths, cell_angles, velocities, kineticEnergy, potentialEnergy, temperature, alchemicalLambda = trj_allatom.read()
    trj_protein.write(coordinates=coordinates[:, protein_atom_indices], time=time, cell_lengths=cell_lengths, cell_angles=cell_angles)  # Ignoring the other fields for now, TODO.

    filenames_protein.append(filenames_allatom[n_files_protein:])
    del trj_allatom, trj_protein

def concatenate_core17(path, top, output_filename):
    """Concatenate tar bzipped XTC files created by Folding@Home Core17.
    
    Parameters
    ----------
    path : str
        Path to directory containing "results-*.tar.bz2".  E.g. a single CLONE directory.
    top : mdtraj.Topology
        Topology for system
    output_filename : str
        Filename of output HDF5 file to generate.
    
    Notes
    -----
    We use HDF5 because it provides an easy way to store the metadata associated
    with which files have already been processed.
    """
    glob_input = os.path.join(path, "results-*.tar.bz2")
    filenames = glob.glob(glob_input)
    filenames = sorted(filenames, key=keynat)
    
    if len(filenames) <= 0:
        return
    
    trj_file = HDF5TrajectoryFile(output_filename, mode='a')
    
    try:
        trj_file._create_earray(where='/', name='processed_filenames',atom=trj_file.tables.StringAtom(1024), shape=(0,))
        trj_file.topology = top.topology
    except trj_file.tables.NodeError:
        pass

    for filename in filenames:
        if six.b(filename) in trj_file._handle.root.processed_filenames:  # On Py3, the pytables list of filenames has type byte (e.g. b"hey"), so we need to deal with this via six.
            print("Already processed %s" % filename)
            continue
        with enter_temp_directory():
            print("Processing %s" % filename)
            archive = tarfile.open(filename, mode='r:bz2')
            archive.extract("positions.xtc")
            trj = md.load("positions.xtc", top=top)

            for frame in trj:
                trj_file.write(coordinates=frame.xyz, cell_lengths=frame.unitcell_lengths, cell_angles=frame.unitcell_angles, time=frame.time)
            
            trj_file._handle.root.processed_filenames.append([filename])

def concatenate_core17_filenames(path, top_filename, output_filename):
    """Concatenate tar bzipped XTC files created by Folding@Home Core17.
    This version accepts only filenames and paths.
    
    Parameters
    ----------
    path : str
        Path to directory containing "results-*.tar.bz2".  E.g. a single CLONE directory.
    top_filename : str
        Filepath to read Topology for system
    output_filename : str
        Filename of output HDF5 file to generate.
    
    Notes
    -----
    We use HDF5 because it provides an easy way to store the metadata associated
    with which files have already been processed.
    """

    print(output_filename)

    # Open topology file.
    top = md.load(top_filename % vars())

    # Glob file paths.
    glob_input = os.path.join(path, "results-*.tar.bz2")
    filenames = glob.glob(glob_input)
    filenames = sorted(filenames, key=keynat)
    
    if len(filenames) <= 0:
        del top
        return
    
    trj_file = HDF5TrajectoryFile(output_filename, mode='a')
    
    try:
        trj_file._create_earray(where='/', name='processed_filenames',atom=trj_file.tables.StringAtom(1024), shape=(0,))
        trj_file.topology = top.topology
    except trj_file.tables.NodeError:
        pass
                
    try:
        for filename in filenames:
            if six.b(filename) in trj_file._handle.root.processed_filenames:  # On Py3, the pytables list of filenames has type byte (e.g. b"hey"), so we need to deal with this via six.
                print("Already processed %s" % filename)
                continue
            with enter_temp_directory():
                print("Processing %s" % filename)
                archive = tarfile.open(filename, mode='r:bz2')
                archive.extract("positions.xtc")
                trj = md.load("positions.xtc", top=top)

                for frame in trj:
                    trj_file.write(coordinates=frame.xyz, cell_lengths=frame.unitcell_lengths, cell_angles=frame.unitcell_angles, time=frame.time)
            
                trj_file._handle.root.processed_filenames.append([filename])

                # Clean up.
                del archive, trj

    except RuntimeError:
        print("Cannot munge RUN%d CLONE%d due to damaged XTC." % (run, clone))
    
    # Clean up.
    del top, trj_file
                
def concatenate_ocore(path, top, output_filename):
    """Concatenate XTC files created by Siegetank OCore.
    
    Parameters
    ----------
    path : str
        Path to stream directory containing frame directories /0, /1, /2
        etc.
    top : mdtraj.Topology
        Topology for system
    output_filename : str
        Filename of output HDF5 file to generate.
    
    Notes
    -----
    We use HDF5 because it provides an easy way to store the metadata associated
    with which files have already been processed.
    """
    sorted_folders = sorted(os.listdir(path), key=lambda value: int(value))
    sorted_folders = [os.path.join(path, folder) for folder in sorted_folders]
    
    if len(sorted_folders) <= 0:
        return
    
    trj_file = HDF5TrajectoryFile(output_filename, mode='a')
    
    try:
        trj_file._create_earray(where='/', name='processed_folders',atom=trj_file.tables.StringAtom(1024), shape=(0,))
        trj_file.topology = top.topology
    except trj_file.tables.NodeError:
        pass
    
    for folder in sorted_folders:
        if six.b(folder) in trj_file._handle.root.processed_folders:  # On Py3, the pytables list of filenames has type byte (e.g. b"hey"), so we need to deal with this via six.
            print("Already processed %s" % folder)
            continue
        print("Processing %s" % folder)
        xtc_filename = os.path.join(folder, "frames.xtc")
        trj = md.load(xtc_filename, top=top)
        
        for frame in trj:
            trj_file.write(coordinates=frame.xyz, cell_lengths=frame.unitcell_lengths, cell_angles=frame.unitcell_angles, time=frame.time)
        
        trj_file._handle.root.processed_folders.append([folder])
            
