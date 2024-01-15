# Masked Comparison Attack
This repository contains data and code used for [redacted].

## Implementation
The implementation can be found under ``recovery/`` and all files require this to be the current working directory.
To run the implementation, first setup a virtual environment (``setup.sh`` will do that), install the requirements in ``requirements.txt``, unpack the traces, and then change to ``recovery/``.
This means, run

```
source setup.sh
cd traces
./unpack_all
cd ../recovery
./main --help
./main physical --help
./main simulation --help
```

or

```
python3 -m venv .env
source .env/bin/activate
pip install -r requirements
cd traces
./unpack_all
cd ../recovery
./main --help
./main physical --help
./main simulation --help
```

The main.py file contains both simulations and physical attack.
To run a simulation use ``./main simulation`` and to run the physical attacks use ``./main physical``.
Attacks can be started with ``--attack [all|template|vertical|horizontal]``.
The switch ``--help`` shows all possible options and can be called on both ``physical`` and ``simulation``.

### Plots
Plots can be obtained using 

```--plot [all|dist-vertical|dist-horizontal|mean|t-test|manual-pois]```

The raw data for all plots is printed to ``../plot_prints/``.

### Physical Attack
Traces and files with bc-values have to be located in ../traces or specified with ``--trace-file`` and ``--bc-file``.
The remaining options can be obtained using ``--help``.

### Simulation
The simulation can be run using

```./main.py simulation --attack all --shares [nshares] --sigma [list of sigmas divided by whitespaces]```.

The remaining options can be obtained using ``--help``.

### Reproducing the Results
To reproduce the results in Table 1, run

```./main.py physical --shares [nshares]```

where ``nshares=3`` or ``nshares=4``.
To reproduce the results in Table 2, run

```./main.py simulation --attack horizontal --ntraces 50 --bcs 5 --shares 4 --sigmas 1 2 5 10 15 20```.

Note that passing ``--ntraces 50`` results in ``2*50`` traces being samples/used of which half are decryption failures and half are decryption successes.

### Manual POI Finding
To manually search for points of interest with potential locations ranging from [start] to [end], use

```./main.py physical --plot manual-pois --select-manual-poi --manual-poi-range [start] [end]```

## Contact
[redacted]
