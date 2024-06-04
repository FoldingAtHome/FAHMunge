"""
Support for core21/core22 projects in ws9 and earlier format.

"""
##############################################################################
# imports
##############################################################################

from __future__ import print_function, division
import os, os.path
import glob
import tarfile
from mdtraj.formats.hdf5 import HDF5TrajectoryFile
import mdtraj as md
import tables
from mdtraj.utils.contextmanagers import enter_temp_directory
from mdtraj.utils import six
from natsort import natsorted
import subprocess
import tempfile
import shutil
import glob
import time
import copy
import sys
import re

################################################################################
# ws9 core21 support
################################################################################

import signal
import time

class SignalHandler:
    """
    """
    terminate = False
    def __init__(self):
        signal.signal(signal.SIGINT, self.exit_gracefully)
        signal.signal(signal.SIGTERM, self.exit_gracefully)

    def exit_gracefully(self, signum, frame):
        self.terminate = True

def list_core21_result_packets(clone_path):
    """
    Create an ordered list of core21 result packets in the specified CLONE path.

    Parameters
    ----------
    clone_path : str
        Path to CLONE directory containing ws8/ws9 WUs

    Returns
    -------
    folders : list of str
        List of core21 result packets, either tarballs (ws7/8) or directories (ws9),
        sorted in order of increasing FRAME number.

    """
    # Get a list of all files in the directory
    file_list = os.listdir(clone_path)

    # Create a list of result packets
    result_packets = list()

    # Iterate over possible result packets
    max_frames = len(file_list)
    for frame_index in range(max_frames):
        # Prioritize uncompressed packets over compressed packets
        result_packet = os.path.join(clone_path, 'results%d' % frame_index)
        if os.path.exists(result_packet):
            result_packets.append(result_packet)
            continue

        # Check if compressed packet exists
        result_packet = os.path.join(clone_path, 'results-%03d.tar.bz2' % frame_index)
        if os.path.exists(result_packet):
            result_packets.append(result_packet)
            continue

        # No more contiguous frames; terminate
        break

    return result_packets

def ensure_result_packet_is_decompressed(result_packet, topology, atom_indices=None, chunksize=10, delete_on_unpack=False, xml_handling='compress'):
    """
    Ensure that the specified result packet is decompressed.

    If this is a ws7/ws8 compressed result packet, safely convert it to uncompressed:
    * decompress it into a temporary directory
    * move it into place
    * verify integrity of files
    * unlink (delete) the old result packet if everything looks OK [OPTIONAL]

    If this is a directory, this function returns immediately.

    .. warning: This will irreversibly delete the compressed work packet, replacing
    it with an uncompressed one.

    Parameters
    ----------
    result_packet : str
        Path to original result packet
    topology : mdtraj.Topology
        Topology to use for verifying integrity of trajectory
    atom_indices : list of int, optional, default=None
        Atom indices to read when verifying integrity of trajectory
        If None, all atoms will be read.
    delete_on_unpack : bool, optional, default=True
        If True, will delete old ws8-style .tar.bz2 files after they have been unpacked.
        WARNING: THIS COULD BE DANGEROUS
    xml_handling : bool, optional, default='compress'
        If 'compress', will compress XML files after unpacking them;
        if 'delete', will delete XML files after unpacking them.
    chunksize : int, optional, default=10
        Number of frames to read each call to mdtraj.iterload for verifying trajectory integrity

    Returns
    -------
    result_packet : str
        Path to new result packet directory

    """
    # Return if this is just a directory
    if os.path.isdir(result_packet):
        return result_packet

    # If this is a tarball, extract salient information.
    # Format: results-002.tar.bz2
    absfilename = os.path.abspath(result_packet)
    (basepath, filename) = os.path.split(absfilename)
    pattern = r'results-(\d+).tar.bz2'
    if not re.match(pattern, filename):
        raise Exception("Compressed results packet filename '%s' does not match expected format (results-001.tar.bz2)" % result_packet)
    frame_number = int(re.match(pattern, filename).group(1))

    # Extract frames from trajectory in a temporary directory
    print("      Extracting %s" % result_packet)
    with enter_temp_directory():
        # Create target directory
        extracted_archive_directory = tempfile.mkdtemp()

        # Extract all contents
        archive = tarfile.open(absfilename, mode='r:bz2')
        archive.extractall(path=extracted_archive_directory)

        # Compress or delete XML files
        xml_filenames = glob.glob('%s/*.xml' % extracted_archive_directory)
        if xml_handling == 'compress':
            for filename in xml_filenames:
                print("      Compressing %s" % os.path.basename(filename))
                subprocess.call(['gzip', filename])
        elif xml_handling == 'delete':
            for filename in xml_filenames:
                print("      Deleting %s" % os.path.basename(filename))
                os.remove(filename)
        else:
            raise Exception('Unknown xml_handling option "%s": Must be "delete" or "compress"' % xml_handling)

        # Create new result packet name
        new_result_packet = os.path.join(basepath, 'results%d' % frame_number)

        # Move directory into place
        shutil.move(extracted_archive_directory, new_result_packet)

        # Verify integrity of archive contents
        xtc_filename = os.path.join(new_result_packet, 'positions.xtc')
        if not os.path.exists(xtc_filename):
            raise Exception("Result packet archive '%s' does not contain positions.xtc; aborting unpacking.")
        try:
            for chunk in md.iterload(xtc_filename, top=topology, atom_indices=atom_indices, chunk=chunksize):
                pass
        except Exception as e:
            msg = "Result packet archive '%s' failed trajectory integrity check; aborting unpacking.\n"
            msg += str(e)
            raise Exception(msg)

        # Cleanup archive object
        del archive

        if delete_on_unpack:
            # Remove archive permanently
            print("      Permanently removing %s" % absfilename)
            os.unlink(absfilename)

        # Return updated result packet directory name
        return new_result_packet

def process_core21_clone(clone_path, topology_filename, processed_trajectory_filename, atom_selection_string, terminate_event=None, delete_on_unpack=False, xml_handling='compress', chunksize=10, signal_handler=None):
    """
    Process core21 result packets in a CLONE, concatenating to a specified trajectory.
    This will append to the specified trajectory if it already exists.

    Note
    ----
    * ws9 stores core21 result packets in an uncompressed directory. Original result packets are left untouched.
    * ws8 and earlier versions store result packets in compressed archives; this method will safely unpack them and remove the original compressed files.
    * An exception will be raised if something goes wrong with processing. The calling process will have to catch this and abort CLONE processing.

    Parameters
    ----------
    clone_path : str
        Source path to CLONE data directory
    topology_filename : str
        Path to PDB or other file containing topology information
    processed_trajectory_filename : str
        Path to concatenated stripped trajectory
    atom_selection_string : str
        MDTraj DSL specifying which atoms should be stripped from source WUs.
    terminate_event : multiprocessing.Event, optional, default=None
        If specified, will terminate early if terminate_event.is_set() is True
    delete_on_unpack : bool, optional, default=True
        If True, will delete old ws8-style .tar.bz2 files after they have been unpacked.
        WARNING: THIS COULD BE DANGEROUS
    xml_handling : bool, optional, default='compress'
        If 'compress', will compress XML files after unpacking them;
        if 'delete', will delete XML files after unpacking them.
    chunksize : int, optional, default=10
        Chunksize (in number of frames) to use for mdtraj.iterload reading of trajectory
    signal_handler : SignalHandler, optional, default=None
        If None, a new SignalHandler object will be created.

    TODO
    ----
    * Add unpacking step to support ws9
    * Include a safer way to substitute vars()

    """
    # Check for early termination since topology reading might take a while
    if terminate_event and terminate_event.is_set():
        return

    MAX_FILEPATH_LENGTH = 1024 # MAXIMUM FILEPATH LENGTH; this may be too short for some installations

    if not signal_handler:
        signal_handler = SignalHandler()

    # Read the topology for the source WU
    # TODO: Only read topology if we have not processed all the WU packets
    # TODO: Use LRU cache to cache work_unit_topology based on filename
    print('Reading topology from %s for clone %s...' % (topology_filename, clone_path))
    top = md.load(topology_filename)
    work_unit_topology = copy.deepcopy(top.topology) # extract topology
    del top # close file

    # Check for early termination since topology reading might take a while
    if terminate_event and terminate_event.is_set():
        return

    # Determine atoms that will be written to trajectory
    atom_indices = work_unit_topology.select(atom_selection_string)

    # Create a new Topology for the atom subset to be written to the trajectory
    trajectory_topology = work_unit_topology.subset(atom_indices)

    # Glob file paths and return result files in sequential order.
    result_packets = list_core21_result_packets(clone_path)

    # Return if there are no WUs to process
    if len(result_packets) <= 0:
        return

    # Open trajectory for appending
    trj_file = HDF5TrajectoryFile(processed_trajectory_filename, mode='a')

    # Initialize new trajectory with topology and list of processed WUs if they are absent
    try:
        # TODO: Switch from pytables StringAtom to arbitrary-length string
        # http://www.pytables.org/usersguide/datatypes.html
        trj_file._create_earray(where='/', name='processed_folders',atom=trj_file.tables.StringAtom(MAX_FILEPATH_LENGTH), shape=(0,))
        trj_file.topology = trajectory_topology # assign topology
    except trj_file.tables.NodeError:
        pass

    # Process each WU, checking whether signal has been received after each.
    for result_packet in result_packets:
        # Stop processing if signal handler indicates we should terminate
        if (signal_handler.terminate) or (terminate_event and terminate_event.is_set()):
            break

        # Skip this WU if we have already processed it
        if six.b(result_packet) in trj_file._handle.root.processed_folders:  # On Py3, the pytables list of filenames has type byte (e.g. b"hey"), so we need to deal with this via six.
            continue

        # If the result packet is compressed, decompress it and return the new directory name
        result_packet = ensure_result_packet_is_decompressed(result_packet, work_unit_topology, delete_on_unpack=delete_on_unpack, xml_handling=xml_handling)

        # Check that we haven't violated our filename length assumption
        if len(result_packet) > MAX_FILEPATH_LENGTH:
            msg = "Filename is longer than hard-coded MAX_FILEPATH_LENGTH limit (%d > %d). Increase MAX_FILEPATH_LENGTH and re-install." % (len(result_packet), MAX_FILEPATH_LENGTH)
            print(msg)
            raise Exception(msg)

        # Stop processing if signal handler indicates we should terminate
        if (signal_handler.terminate) or (terminate_event and terminate_event.is_set()):
            break

        # Process the work unit
        # TODO: Write to logger instead of printing to terminal
        # TODO: We could conceivably also check for early termination in the chunk loop if we carefully track the last chunk processed as well.
        print("   Processing %s" % result_packet)
        xtc_filename = os.path.join(result_packet, "positions.xtc")
        for chunk in md.iterload(xtc_filename, top=work_unit_topology, atom_indices=atom_indices, chunk=chunksize):
            trj_file.write(coordinates=chunk.xyz, cell_lengths=chunk.unitcell_lengths, cell_angles=chunk.unitcell_angles, time=chunk.time)
        # Record that we've processed the WU
        trj_file._handle.root.processed_folders.append([result_packet])

    # Sync the trajectory file to flush all data to disk
    trj_file.close()

    # Make sure we tell everyone to terminate if we are terminating
    if signal_handler.terminate and terminate_event:
        terminate_event.set()
