#!/usr/bin/env python3

# ##### DESCRIPTION ######
# By default this file is to be called as cmdline application
# However, you may work with it by hand by editing the global variables
# and switching out the main function at the very end of the file
# to main_simulation or main_physical
##########################

import sys
import argparse
import numpy as np
from datetime import datetime
from poi import PoisCollection
from util import read_traces, find_poi_manual
from attack import attack, attack_one_trace
from plotting import plot_diff_means, plot_t_test_fail_nfail, plot_distribution_bc_horizontal, plot_distribution_bc_vertical
from simulation import Simulator

# #### SHARED PARAMS ####
NSHARES = 4
NTRACES = 500
NTRACES_PROFILE = 500
NUM_BCS = 1
POIS_PER_BIT = 1

PERFORM_PLOT_HORIZONTAL = False
PERFORM_PLOT_VERTICAL = False

# #### PHYSICAL PARAMS ####
DECIMATE = 1
OPT_LEVEL = 2

PERFORM_TEMPLATE_ONE_TRACE = False
PERFORM_TEMPLATE = False

PERFORM_VERTICAL_ONE_TRACE = False
PERFORM_VERTICAL = False

PERFORM_HORIZONTAL_ONE_TRACE = False
PERFORM_HORIZONTAL = False

PERFORM_PLOT_MEAN = False
PERFORM_T_TEST = False
PERFORM_PLOT_PLOI_FINDING = False

SEPARATE_TEMPLATE = False

TEST_MODE = True

TRACE_FILE = None
BC_FILE = None

TRACE_INDEX = 0
MANUAL_POI_RANGE = None
NO_POI_FINDING = False

# #### SIMULATION PARAMS ####
SIGMAS = [4.0]
SEED = 42
EVAL_TRACES = None


def main():
    shared_parser = argparse.ArgumentParser(add_help=False)
    shared_parser.add_argument('--attack', choices=['all', 'template', 'vertical', 'horizontal'], nargs="+", default=[])
    shared_parser.add_argument('--plot', choices=['all', 'mean', 't-test', 'manual-pois', 'dist-horizontal', 'dist-vertical'], nargs="+", default=[])
    shared_parser.add_argument('--bcs', type=int, default=1)
    shared_parser.add_argument('--pois', choices=range(1, 4), type=int, default=1)
    shared_parser.add_argument('--shares', choices=range(2, 20), type=int, default=4)
    shared_parser.add_argument('--ntraces', type=int, default=500)
    shared_parser.add_argument('--ntraces-profile', type=int, default=500)

    parser = argparse.ArgumentParser(prog='Attacking Masked Comparisons')
    subparsers = parser.add_subparsers(required=True, dest='simulation_or_physical')

    parser_physical = subparsers.add_parser('physical', parents=[shared_parser])
    parser_physical.add_argument('--attack-one-trace', choices=['all', 'template', 'vertical', 'horizontal'], nargs="+", default=[])
    parser_physical.add_argument('--decimate', choices=[1, 10], type=int, default=1)
    parser_physical.add_argument('--optlevel', choices=[2, 3], type=int, default=2)
    parser_physical.add_argument('--separate-template', action="store_true")
    parser_physical.add_argument('--no-test', action="store_true")
    parser_physical.add_argument('--trace-file', type=str, default=None)
    parser_physical.add_argument('--bc-file', type=str, default=None)
    parser_physical.add_argument('--trace-index', type=int, default=0)
    parser_physical.add_argument('--manual-poi-range', type=int, default=[450, 450], nargs=2)
    parser_physical.add_argument('--select-manual-poi', action='store_true')

    parser_simulation = subparsers.add_parser('simulation', parents=[shared_parser])
    parser_simulation.add_argument("--seed", type=int, help="Simulation seed", default=42)
    parser_simulation.add_argument("-s", "--sigmas", type=float, help="Noise level in standard deviation sigma", default=[5.0], nargs="+")
    parser_simulation.add_argument("-e", "--eval-traces", type=int, help="Number of traces used for evaluation", default=None)

    args = parser.parse_args()
    args_dict = vars(args)

    global SIGMAS
    SIGMAS = args_dict.get('sigmas')
    global EVAL_TRACES
    EVAL_TRACES = args_dict.get('eval_traces')
    global SEED
    SEED = args_dict.get('seed')

    global TRACE_FILE
    TRACE_FILE = args_dict.get('trace_file')
    global BC_FILE
    BC_FILE = args_dict.get('trace_file')

    global TRACE_INDEX
    TRACE_INDEX = args_dict.get('trace_index')

    global NTRACES
    NTRACES = args_dict.get('ntraces')
    global NTRACES_PROFILE
    NTRACES_PROFILE = args_dict.get('ntraces_profile')
    global OPT_LEVEL
    OPT_LEVEL = args_dict.get('optlevel')
    global DECIMATE
    DECIMATE = args_dict.get('decimate')
    global NSHARES
    NSHARES = args_dict.get('shares')
    global POIS_PER_BIT
    POIS_PER_BIT = args_dict.get('pois')
    global NUM_BCS
    NUM_BCS = args_dict.get('bcs')

    if NTRACES_PROFILE % 2 != 0:
        print("Error: Number of profile traces have to be even for implementation reasons.")
        exit(1)

    global SEPARATE_TEMPLATE
    SEPARATE_TEMPLATE = args_dict.get('separate_template')
    global TEST_MODE
    TEST_MODE = not args_dict.get('no_test')
    global MANUAL_POI_RANGE
    MANUAL_POI_RANGE = args_dict.get('manual_poi_range')
    global NO_POI_FINDING
    NO_POI_FINDING = args_dict.get('select_manual_poi')

    plot_list = args_dict.get('plot')
    if plot_list is None:
        plot_list = []
    plot_all = 'all' in plot_list
    if 'mean' in plot_list or plot_all:
        global PERFORM_PLOT_MEAN
        PERFORM_PLOT_MEAN = True
    if 't-test' in plot_list or plot_all:
        global PERFORM_T_TEST
        PERFORM_T_TEST = True
    if 'manual-pois' in plot_list or plot_all or NO_POI_FINDING:
        global PERFORM_PLOT_PLOI_FINDING
        PERFORM_PLOT_PLOI_FINDING = True
        print("WARNING: THIS OPTION REQUIRES ADJUSTING THE PARAMETERS BY HAND")
    if 'dist-horizontal' in plot_list or plot_all:
        global PERFORM_PLOT_HORIZONTAL
        PERFORM_PLOT_HORIZONTAL = True
    if 'dist-vertical' in plot_list or plot_all:
        global PERFORM_PLOT_VERTICAL
        PERFORM_PLOT_VERTICAL = True

    attack_list = args_dict.get('attack')
    if attack_list is None:
        attack_list = []
    perform_all = 'all' in attack_list
    if 'template' in attack_list or perform_all:
        global PERFORM_TEMPLATE
        PERFORM_TEMPLATE = True
    if 'vertical' in attack_list or perform_all:
        global PERFORM_VERTICAL
        PERFORM_VERTICAL = True
    if 'horizontal' in attack_list or perform_all:
        global PERFORM_HORIZONTAL
        PERFORM_HORIZONTAL = True

    attack_one_trace_list = args_dict.get('attack_one_trace')
    if attack_one_trace_list is None:
        attack_one_trace_list = []
    perform_all_one_trace = 'all' in attack_one_trace_list
    if 'template' in attack_one_trace_list or perform_all_one_trace:
        global PERFORM_TEMPLATE_ONE_TRACE
        PERFORM_TEMPLATE_ONE_TRACE = True
    if 'vertical' in attack_one_trace_list or perform_all_one_trace:
        global PERFORM_VERTICAL_ONE_TRACE
        PERFORM_VERTICAL_ONE_TRACE = True
    if 'horizontal' in attack_one_trace_list or perform_all_one_trace:
        global PERFORM_HORIZONTAL_ONE_TRACE
        PERFORM_HORIZONTAL_ONE_TRACE = True

    print(f"Received arguments: {args_dict}")
    if args.simulation_or_physical == 'physical':
        main_physical()
    elif args.simulation_or_physical == 'simulation':
        main_simulation()
    else:
        raise ValueError


def print_base_settings(sim):
    print()
    if sim:
        print("#"*10 + "#"*18 + "#"*10)
        print("#"*10 + " SIMULATED ATTACK " + "#"*10)
        print("#"*10 + "#"*18 + "#"*10)
    else:
        print("#"*10 + "#"*17 + "#"*10)
        print("#"*10 + " PHYSICAL ATTACK " + "#"*10)
        print("#"*10 + "#"*17 + "#"*10)
    print()
    print("#"*10 + " SETTINGS " + "#"*10)
    print(f"{NSHARES=}, {DECIMATE=}, {OPT_LEVEL=}, {NTRACES=}, {NTRACES_PROFILE=}, {POIS_PER_BIT=}")
    print("#"*10 + "#"*10 + "#"*10)
    print()


def read_all_traces(directory="traces"):
    print("Reading traces..")
    if TRACE_FILE is None:
        file = f"../{directory}/traces-order-{NSHARES}-decimate-{DECIMATE}-O{OPT_LEVEL}-{NTRACES}.bin"
    else:
        file = TRACE_FILE
    if BC_FILE is None:
        file_bc = f"../{directory}/bcs-order-{NSHARES}.bin"
    else:
        file_bc = BC_FILE
    print(f"Trace file: {file}")
    print(f"BC file: {file_bc}")
    try:
        traces, bcs = read_traces(NTRACES, NSHARES, file=file, file_bc=file_bc, ignore_bc=SEPARATE_TEMPLATE)
        if SEPARATE_TEMPLATE:
            file_profile = f"../{directory}/traces-order-{NSHARES}-decimate-{DECIMATE}-O{OPT_LEVEL}-{NTRACES_PROFILE}.bin"
            file_bc_profile = file_bc
            print(f"Profile trace file: {file_profile}")
            print(f"Profile BC file: {file_bc_profile}")
            traces_profile, bcs_profile = read_traces(NTRACES_PROFILE, NSHARES, file=file_profile, file_bc=file_bc_profile, ignore_bc=False)
            print("Attack traces for profiled and non-profiled attacks are the same.")
            print(f"Profiled attacks use additional traces from {file_profile}.")
            traces_attack = traces
            bcs_attack = bcs
        else:
            print(f"First {NTRACES_PROFILE} will be used for profiling (i.e., for finding POIs and the template attacks).")
            print(f"The remaining {2*NTRACES-NTRACES_PROFILE} traces will be used for the template attacks.")
            print("All traces will be used in the vertical and horizontal attacks.")
            traces_profile = traces[:NTRACES_PROFILE]
            bcs_profile = bcs[:NTRACES_PROFILE]
            traces_attack = traces[NTRACES_PROFILE:]
            bcs_attack = bcs[NTRACES_PROFILE:]
    except OSError as e:
        print(f"Unable to read trace or bc file: {e}", file=sys.stderr)
        print("Are you sure the traces for the selected setting exist?")
        return
    print()
    return traces_profile, bcs_profile, traces_attack, bcs_attack, traces, bcs


def perform_plots(traces, bcs, pois, trace_idx):
    if PERFORM_PLOT_PLOI_FINDING:
        print("Plotting manual poi finding..")
        print("WARNING: THIS OPTION REQUIRES ADJUSTING THE PARAMETERS BY HAND")
        if pois is not None:
            print("POIs: " + str(list(map(lambda x: x.trace_locs[0], pois.get_pois_list()))))
        if NO_POI_FINDING:
            pois = find_poi_manual(traces, start=MANUAL_POI_RANGE[0], end=MANUAL_POI_RANGE[1], sel=True, shares=NSHARES)
        else:
            find_poi_manual(traces, start=MANUAL_POI_RANGE[0], end=MANUAL_POI_RANGE[1])
        print()

    if PERFORM_PLOT_MEAN:
        print("Plotting..")
        plot_diff_means(traces, bcs, pois=pois)
        print()

    if PERFORM_T_TEST:
        print("Plotting t-test..")
        plot_t_test_fail_nfail(traces, pois)
        print()

    ###########
    if PERFORM_PLOT_HORIZONTAL:
        print("Plotting horizontal..")
        if bcs is not None:
            plot_distribution_bc_horizontal(traces[trace_idx], pois, bcs[trace_idx])
        else:
            plot_distribution_bc_horizontal(traces[trace_idx], pois)
        print()
    ###########

    ###########
    if PERFORM_PLOT_VERTICAL:
        print("Plotting vertical..")
        plot_distribution_bc_vertical(traces, pois, bcs)
        print()
    ###########
    return pois


def perform_attacks(pois, traces, traces_attack, bcs_attack, traces_profile, bcs_profile, trace_idx):
    ###########
    if PERFORM_TEMPLATE_ONE_TRACE:
        def templ_func_0():
            return pois.compute_template_from_bc(traces_profile, bcs_profile)
        attack_one_trace(traces_attack, bcs_attack, trace_idx, pois, name=f"Single trial template {NSHARES}-{NTRACES}-{DECIMATE}", template_function=templ_func_0)
    ###########

    ###########
    if PERFORM_TEMPLATE:
        def templ_func_1():
            return pois.compute_template_from_bc(traces_profile, bcs_profile)
        attack(traces_attack, pois, NSHARES, f"Template attack {NSHARES}-{NTRACES}-{DECIMATE}", test_mode=TEST_MODE, template_function=templ_func_1)
    ###########

    ###########
    if PERFORM_VERTICAL_ONE_TRACE:
        def templ_func_2():
            return pois.compute_vertical_auto_template(traces)
        attack_one_trace(traces_attack, bcs_attack, trace_idx, pois, name=f"Single trial vertical {NSHARES}-{NTRACES}-{DECIMATE}", template_function=templ_func_2)
    ###########

    ###########
    if PERFORM_VERTICAL:
        def templ_func_3():
            return pois.compute_vertical_auto_template(traces)
        attack(traces_attack, pois, NSHARES, f"Vertical attack: {NSHARES}-{NTRACES}-{DECIMATE}", test_mode=TEST_MODE, template_function=templ_func_3)
    ###########

    ###########
    if PERFORM_HORIZONTAL_ONE_TRACE:
        def templ_func_4():
            return pois.compute_horizontal_auto_template(traces[trace_idx])
        attack_one_trace(traces_attack, bcs_attack, trace_idx, pois, name=f"Single trial horizontal {NSHARES}-{NTRACES}-{DECIMATE}", template_function=templ_func_4)
    ###########

    ###########
    if PERFORM_HORIZONTAL:
        attack(traces, pois, NSHARES, f"Horizontal attack: {NSHARES}-{NTRACES}-{DECIMATE}", build_single_trace_template=PoisCollection.compute_horizontal_auto_template, test_mode=TEST_MODE, template_function=None)
    ###########


def main_physical():
    print_base_settings(sim=False)

    trace_idx = TRACE_INDEX

    traces_profile, bcs_profile, traces_attack, bcs_attack, traces, bcs = read_all_traces()

    pois = None
    ###########
    if not NO_POI_FINDING:
        print("Finding pois..")
        pois = PoisCollection.find_all_pois(traces_profile, bcs_profile, num_pois_per_bit=POIS_PER_BIT)
        print(f"Found {pois.get_num_pois_per_bit()} POIs per bit for {pois.get_num_bcs()} BCs and {pois.get_num_shares()} shares.")
        print()
    ###########

    pois_2 = perform_plots(traces, bcs, pois, trace_idx)
    if NO_POI_FINDING:
        pois = pois_2
    perform_attacks(pois, traces, traces_attack, bcs_attack, traces_profile, bcs_profile, trace_idx)


def main_simulation():
    results = {}
    now = datetime.now()
    dt_string = now.strftime("%d%m%Y-%h%m%s")

    print_base_settings(sim=True)

    if POIS_PER_BIT > 1:
        print("ERROR: MULTIPLE POIs ARE NOT MODELED CORRECTLY. ABORTING.")
        exit(1)
    for sigma in SIGMAS:
        print(f"Seeding with {SEED}.")
        np.random.seed(SEED)
        print(f"Creating simulator with {sigma=}, {NUM_BCS=}, {NSHARES=}, {POIS_PER_BIT=}, {2*NTRACES=}, {EVAL_TRACES=}..")
        sim = Simulator(sigma, NUM_BCS, NSHARES, POIS_PER_BIT)
        print("Finding POIs..")
        sim.find_pois()
        print("Recording traces..")
        sim.record_traces(2*NTRACES)
        sim.finish_recording_phase()
        if PERFORM_PLOT_VERTICAL:
            print("Plotting vertical..")
            plot_distribution_bc_vertical(sim.traces, sim.pois, sim.bcs)
        if PERFORM_PLOT_HORIZONTAL:
            print("Plotting horizontal..")
            plot_distribution_bc_horizontal(sim.traces[0], sim.pois, sim.bcs[0])
        print("Executing attacks..")
        res = sim.execute_attacks(number_of_traces=EVAL_TRACES, profile_traces=NTRACES_PROFILE, template=PERFORM_TEMPLATE, vertical=PERFORM_VERTICAL, horizontal=PERFORM_HORIZONTAL)
        results[sigma] = res
    print(results)
    results_file = f"../results/results_{dt_string}.txt"
    try:
        with open(results_file, 'w') as f:
            f.write(str(results))
    except OSError as e:
        print(f"Unable to write {results_file}: {e}", file=sys.stderr)
        print("Are you sure the results directory exists?")

    return results


if __name__ == "__main__":
    main()
