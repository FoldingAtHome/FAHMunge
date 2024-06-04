[//]: # (Badges)
[![Build Status](https://travis-ci.org/choderalab/fahmunge.png)](https://travis-ci.org/choderalab/fahmunge)
[![codecov](https://codecov.io/gh/choderalab/fahmunge/branch/master/graph/badge.svg)](https://codecov.io/gh/choderalab/fahmunge/branch/master)
[![Anaconda Cloud Badge](https://anaconda.org/omnia/fahmunge/badges/version.svg)](https://anaconda.org/omnia/fahmunge)
[![Downloads](https://anaconda.org/omnia/fahmunge/badges/downloads.svg)](https://anaconda.org/omnia/fahmunge)

## fahmunge

A tool for automated processing of Folding@home data to produce [mdtraj](http://mdtraj.org)-compatible trajectory sets.

#### Authors
* Kyle A. Beauchamp
* John D. Chodera
* Steven Albanese
* Rafal Wiewiora

#### Supported FAH cores
* 0x17 [OpenMM]
* 0x18 [OpenMM]
* 0x21 [OpenMM 6.3]

#### Installation

The easiest way to install `fahmunge` and its dependencies is via `conda` (preferably [`miniconda`](http://conda.pydata.org/miniconda.html)):
```bash
conda install --yes -c omnia fahmunge
```

#### Usage

##### Basic Usage

Basic usage simply specifies a project CSV file and an output path for the munged data:
```bash
munge-fah-data --projects projects.csv --outpath /data/choderalab/fah/munged3 --nprocesses 16 --validate
```
The metadata for FAH is a CSV file located here on `choderalab` FAH servers:
```
/data/choderalab/fah/Software/FAHMunge/projects.csv
```
This file specifies the project number, the location of the FAH data, a reference PDB file (or files) to be used for munging, and the MDTraj DSL topology selection to be used for extracting solute coordinates of interest.

For example:
```
project,location,pdb,topology_selection
"10491","/home/server.140.163.4.245/server2/data/SVR2359493877/PROJ10491/","/home/server.140.163.4.245/server2/projects/GPU/p10491/topol-renumbered-explicit.pdb","not water"
"10492","/home/server.140.163.4.245/server2/data/SVR2359493877/PROJ10492/","/home/server.140.163.4.245/server2/projects/GPU/p10492/topol-renumbered-explicit.pdb","not (water or resname NA or resname CL)"
"10495","/home/server.140.163.4.245/server2/data/SVR2359493877/PROJ10492/","/home/server.140.163.4.245/server2/projects/GPU/p10495/MTOR_HUMAN_D0/RUN%(run)d/system.pdb","not (water or resname NA or resname CL)"
```
`pdb` points the pipeline toward a PDB file to look at for numbering atoms in the munged data.
The top two lines are examples of using a single PDB for all RUNs in the project.
The third line shows how to use a different PDB for each RUN.
`%(run)d` is substituted by the run number via `filename % vars()` in Python, which allows run numbers or other local Python variables to be substituted.
Substitution is only performed on a per-run basis, not per-clone.

The projects CSV file will undergo minimal validation automatically to make sure all data and file paths can be found.

##### Advanced Usage

More advanced usage allows additional arguments to be specified:
* `--nprocesses <NPROCESSES>` will parallelize munging by RUN using `multiprocessing` if `NPROCESSES > 1` is specified.  By default, `NPROCESSES = 1`.
* `--time <TIME_LIMIT>` specifies that munging should move on to another phase or project after the given time limit (in seconds) is reached, once it is safe to move on.  This is useful for ensuring that some munging occurs on all projects of interest every day.
* `--verbose` will produce verbose output
* `--maxits <MAXITS>` will cause the munging pipeline to run for the specified number of iterations and then exit. This can be useful for debugging. Without specifying this option, munging will run indefinitely.
* `--sleeptime <SLEEPTIME>` will cause munging to sleep for the specified number of seconds if no work was done in this iteration (default:3600).
* `--validate` will validate the choice of `topology_selection` MDTraj DSL topology selection queries to make sure they are valid; note that this may take a significant amount of time, so is optional behavior
* `--xml compress` will compress `.xml` files after unpacking them from old-WS-style result packages to save space; `--xml delete` will delete `.xml` files after unpacking
* `--delete-tarballs` will delete the old `.tar.bz2` WUs after unpacking them from old-WS-style result packages to save space

#### Usage on `choderalab` Folding@home servers

1.  Login to work server using the usual FAH login
2.  Check if script is running (`screen -r -d`).  If True, stop here.
3.  Start a screen session
4.  Run with: `munge-fah-data --projects /data/choderalab/fah/projects.csv --outpath /data/choderalab/fah/munged-data --time 600 --nprocesses 16 --delete-on-unpack --xml compress`
5.  To stop, control c when the script is in the "sleep" phase

#### How it works

Overall Pipeline (Core17/18):

1.  Extract XTC data from `bzip`s
2.  Append all-atom coordinates and filenames to HDF5 file
3.  Extract protein coordinates and filenames from the all-atom HDF5 file into a second HDF5 file

#### Efficiency considerations

The rate limiting step appears to be `bunzip`.  
If we can avoid having the trajectories be double-`bzip`ped by the client, this will speed up things immensely.

#### Nightly syncing to `hal.cbio.mskcc.org`

Munged `no-solvent` data is `rsync`ed nightly from `plfah1` and `plfah2` to `hal.cbio.mskcc.org` via the `choderalab` robot user account to:
```
/cbio/jclab/projects/fah/fah-data/munged
```
This is done via a `crontab`:
```
# kill any rsyncs already in progress
42 00 * * * skill rsync
# munged3
04 01 * * * rsync -av --append-verify --bwlimit=1000 --chmod=g-w,g+r,o-w,o+r server@plfah1.mskcc.org:/data/choderalab/fah/munged2/no-solvent /cbio/jclab/projects/fah/fah-data/munged3 >> $HOME/plfah1-rsync3-no-solvent.log 2>&1
38 02 * * * rsync -av --append-verify --bwlimit=1000 --chmod=g-w,g+r,o-w,o+r server@plfah2.mskcc.org:/data/choderalab/fah/munged2/no-solvent /cbio/jclab/projects/fah/fah-data/munged3 >> $HOME/plfah2-rsync3-no-solvent.log 2>&1
34 03 * * * rsync -av --append-verify --bwlimit=1000 --chmod=g-w,g+r,o-w,o+r server@plfah1.mskcc.org:/data/choderalab/fah/munged2/all-atoms /cbio/jclab/projects/fah/fah-data/munged3 >> $HOME/plfah1-rsync3-all-atoms.log 2>&1
50 03 * * * rsync -av --append-verify --bwlimit=1000 --chmod=g-w,g+r,o-w,o+r server@plfah2.mskcc.org:/data/choderalab/fah/munged2/all-atoms /cbio/jclab/projects/fah/fah-data/munged3 >> $HOME/plfah2-rsync3-all-atoms.log 2>&1
```
To install this `crontab` as the `choderalab` user:
```bash
crontab ~/crontab
```
To list the active `crontab`:
```bash
crontab -l
```
Transfers are logged in the `choderalab` account:
```
plfah1-rsync3-all-atoms.log
plfah1-rsync3-no-solvent.log
plfah2-rsync3-all-atoms.log
plfah2-rsync3-no-solvent.log
```

#### Acknowledgements

Project skeleton based on the
[Computational Chemistry Python Cookiecutter](https://github.com/choderalab/cookiecutter-python-comp-chem)
