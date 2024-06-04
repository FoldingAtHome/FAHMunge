import time
import os
import glob
import argparse
import sys
import collections
import datetime
import pandas as pd
import mdtraj as md

import fahmunge

# Reads in a list of project details from a CSV file with Core17/18 FAH projects and munges them.

def setup_worker(terminate_event, delete_on_unpack, xml_handling):
    global global_terminate_event
    global_terminate_event = terminate_event
    global global_delete_on_unpack
    global_delete_on_unpack = delete_on_unpack
    global xml_handling
    global_xml_handling = xml_handling

def worker(args):
    return fahmunge.core21.process_core21_clone(*args, terminate_event=global_terminate_event, delete_on_unpack=global_delete_on_unpack, compress_xml=global_compress_xml)

def main():
    description = 'Munge FAH data'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-p', '--projects', metavar='PROJECTFILE', dest='projectfile', action='store', type=str, default=None,
        help='CSV file containing (project,filepath,pdbfile) tuples')
    parser.add_argument('--validate', dest='validate_topology_selection', action='store_true', default=False,
        help='Validate topology_selection in CSV projects file is valid')
    parser.add_argument('-o', '--outpath', metavar='OUTPATH', dest='output_path', action='store', type=str, default=None,
        help='Output pathname for munged data')
    parser.add_argument('-n', '--nprocesses', metavar='NPROCESSES', dest='nprocesses', action='store', type=int, default=1,
        help='Number of threads to use (default: 1)')
    parser.add_argument('-d', '--debug', dest='debug', action='store_true', default=False,
        help='Run in serial mode and turn on debug output')
    parser.add_argument('-z', '--delete-tarballs', dest='delete_on_unpack', action='store_true', default=False,
        help='Delete original results-###.tar.bz2 after unpacking; WARNING: THIS IS DANGEROUS AND COULD DELETE YOUR PRIMARY DATA.')
    parser.add_argument('-t', '--time', metavar='TIME', dest='time_limit', action='store', type=int, default=None,
        help='Process each project for no more than specified time (in seconds) before moving on to next project')
    parser.add_argument('-m', '--maxits', metavar='MAXITS', dest='maximum_iterations', action='store', type=int, default=None,
        help='Perform specified number of iterations and exist (default: no limit, process indefinitely)')
    parser.add_argument('-s', '--sleeptime', metavar='SLEEPTIME', dest='sleep_time', action='store', type=int, default=0,
        help='Sleep for specified time (in seconds) between iterations (default: 0)')
    parser.add_argument('-v', '--version', action='store_true', default=False,
        help='Print version information and exit')
    parser.add_argument('-x', '--xml', dest='xml_handling', action='store', default='compress',
        help='How to handle old WS style .xml files: "compress" (default) or "delete"')
    args = parser.parse_args()

    if args.version:
        print(fahmunge.__version__)
        sys.exit(0)

    # Check arguments
    if args.projectfile is None:
        print('ERROR: projectfile must be specified\n\n')
        parser.print_help()
        sys.exit(1)
    if args.output_path is None:
        print('ERROR: outpath must be specified\n\n')
        parser.print_help()
        sys.exit(1)
    if args.nprocesses <= 0:
        print('ERROR: nprocesses must be positive\n\n')
        parser.print_help()
        sys.exit(1)

    # Read project tuples
    projects = pd.read_csv(args.projectfile, index_col=0)

    # Check that all locations and PDB files exist, raising an exception if they do not (indicating misconfiguration)
    # TODO: Parallelize validation by entry?
    print('Validating contents of project CSV file...')
    for (project, location, pdb, topology_selection) in projects.itertuples():
        # Check project path exists.
        if not os.path.exists(location):
            raise Exception("Project %s: Cannot find data path '%s'. Check that you specified the correct location." % (project, location))
        # Check PDB file(s) exist
        # TODO: Check atom counts match?
        # TODO: Generalize with a generator for iterating over all RUN/CLONEs?
        n_runs, n_clones = fahmunge.automation.get_num_runs_clones(location)
        print("Project %s: %d RUNs %d CLONEs found; topology_selection = '%s'" % (project, n_runs, n_clones, topology_selection))
        if '%' in pdb:
            # perform filename substitution on all RUNs
            pdb_filenames_to_check = list()
            for run in range(n_runs):
                pdb_filename = pdb % vars()
                pdb_filenames_to_check.append(pdb_filename)
        else:
            pdb_filenames_to_check = [ pdb ] # just one filename
        for pdb_filename in pdb_filenames_to_check:
            if not os.path.exists(pdb_filename):
                raise Exception("Project %s: PDB filename specified as '%s' but '%s' was not found. Check that you specified the correct path and PDB files are present." % (project, pdb, pdb_filename))
            if args.validate_topology_selection:
                # Check toplogy selection is valid.
                # TODO: Report on original and stripped atom numbers
                traj = md.load(pdb_filename)
                original_topology = traj.top
                indices = original_topology.select(topology_selection)
                print("  %s : %d atoms; selection '%s' has %d atoms" % (pdb_filename, original_topology.n_atoms, topology_selection, len(indices)))
                if len(indices)==0:
                    raise Exception("topology_selection '%s' matches zero atoms!" % topology_selection)
                del traj, original_topology, indices
    print('All specified paths and PDB files found.')
    print('')

    # Report any special processing requests
    if args.maximum_iterations:
        print('Processing for a total of %d iterations' % args.maximum_iterations)
    if args.time_limit:
        print('Will run for %s seconds and terminate' % args.time_limit)
    if args.sleep_time:
        print('Will sleep for %s seconds between iterations' % args.sleep_time)
    print('')

    # Set signal handling
    signal_handler = fahmunge.core21.SignalHandler()

    # Main processing loop
    iteration = 0
    terminate = False # if True, terminate
    initial_time = time.time()
    while(not terminate):
        # Assemble list of CLONEs to process
        print('----------' * 8)
        print('Iteration %8d : Assembling list of CLONEs to process...' % iteration)
        print(datetime.datetime.now().isoformat())
        print('----------' * 8)
        clones_to_process = collections.deque()
        for (project, project_path, topology_filename, topology_selection) in projects.itertuples():

            print('Project %s' % project)
            print("  location: '%s'" % project_path)
            print("  reference topology file: '%s'" % topology_filename)
            print("  topology selection: '%s'" % topology_selection)

            # Form output path
            output_path = os.path.join(args.output_path, "%s/" % project)

            # Make sure output path exists
            fahmunge.automation.make_path(output_path)

            # Determine number of RUNs and CLONEs
            n_runs, n_clones = fahmunge.automation.get_num_runs_clones(project_path)

            # Compile CLONEs to process
            for run in range(n_runs):
                for clone in range(n_clones):
                    # Get clone source and destination paths
                    clone_path = os.path.join(project_path, "RUN%d" % run, "CLONE%d" % clone)
                    processed_clone_filename = os.path.join(output_path, "run%d-clone%d.h5" % (run, clone))
                    # Form work packet
                    work_args = (clone_path, topology_filename % vars(), processed_clone_filename, topology_selection)
                    # Append work packet
                    clones_to_process.append(work_args)

            # Terminate if instructed
            if signal_handler.terminate:
                print('Signal caught; terminating.')
                exit(1)

        print('There are %d CLONEs to process' % len(clones_to_process))
        print('----------' * 8)
        print('')

        # Munge data in parallel
        print('----------' * 8)
        print('Iteration %8d : Processing %d CLONEs...' % (iteration, len(clones_to_process)))
        print(datetime.datetime.now().isoformat())

        if args.debug:
            print('Using serial debug mode')
            print('----------' * 8)
            for packed_args in clones_to_process:
                fahmunge.core21.process_core21_clone(*packed_args, delete_on_unpack=args.delete_on_unpack, xml_handling=args.xml_handling, signal_handler=signal_handler)
                # Terminate if instructed
                if signal_handler.terminate:
                    print('Signal caught; terminating.')
                    exit(1)
        else:
            # Settings for thread processing
            print('Using %d threads' % args.nprocesses)
            print('----------' * 8)
            from multiprocessing import Pool, Event
            print("Creating thread pool of %d threads..." % args.nprocesses)
            terminate_event = Event()
            pool = Pool(args.nprocesses, setup_worker, (terminate_event, args.delete_on_unpack, args.xml_handling))


            try:
                print("Starting asynchronous map operations...")
                job = pool.map_async(worker, clones_to_process, chunksize=1)

                sleep_interval = 5 # seconds between polling of multiprocessing pool
                while( (not job.ready()) and (not terminate_event.is_set()) ):
                    time.sleep(sleep_interval)
                    # Terminate if maximum time has elapsed.
                    elapsed_time = time.time() - initial_time
                    if args.time_limit and (elapsed_time > args.time_limit):
                        print('Elapsed time (%.1f s) exceeds timeout (%.1f s); signaling jobs to terminate.' % (elapsed_time, args.time_limit))
                        terminate_event.set()
                        terminate = True
                    # Terminate if a signal has been caught
                    if signal_handler.terminate:
                        print('Signal caught; terminating.')
                        terminate_event.set()
                        terminate = True

            except KeyboardInterrupt:
                print("Caught KeyboardInterrupt, safely terminating workers. This may take several minutes. Please be patient to avoid data corruption.")
                # Signal termination
                terminate_event.set()
                terminate = True
                # Close down the multiprocessing pool
                pool.close()
                pool.join()

            except Exception as e:
                print('An exception occurred; terminating...')
                # An exception occurred; terminate.
                print(e)
                raise e

            finally:
                print("Cleaning up...")
                pool.close()
                pool.join()

        # Report completion of iteration
        print('Finished iteration %d.' % iteration)

        # Increment iteration counter
        iteration += 1

        # If time limit has elapsed, terminate
        elapsed_time = time.time() - initial_time
        if args.time_limit and (elapsed_time > args.time_limit):
            print('Elapsed time (%.1f s) exceeds timeout (%.1f s); signaling jobs to terminate.' % (elapsed_time, args.time_limit))
            terminate = True

        # Exit now if specified number of iterations is reached
        if args.maximum_iterations and (iteration >= args.maximum_iterations):
            print('Maximum number of iterations (%d) reached.' % args.maximum_iterations)
            terminate = True

        if terminate:
            return

        # Sleep
        print("Sleeping for %d seconds." % (args.sleep_time))
        time.sleep(args.sleep_time)

        # If time limit has elapsed, terminate
        elapsed_time = time.time() - initial_time
        if args.time_limit and (elapsed_time > args.time_limit):
            print('Elapsed time (%.1f s) exceeds timeout (%.1f s); signaling jobs to terminate.' % (elapsed_time, args.time_limit))
            terminate = True

        # End of iteration
        print('----------' * 8)
        print('')
