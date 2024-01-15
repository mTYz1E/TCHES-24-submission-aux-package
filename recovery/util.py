import sys
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

SIMPLECOMPBITS = 272


def share_value(bc, nshares):
    shares = [np.random.randint(0, 1 << 32) for _ in range(nshares-1)]
    bc_shared = bc
    for share in shares:
        bc_shared ^= share
    shares.append(bc_shared)
    unshared = 0
    for share in shares:
        unshared ^= share
    assert bc == unshared
    return shares


def bits_2(num):
    return [int(b) for b in bin(num)[2:]]


def bits(num):
    bits = []
    for _ in range(32):
        bits.append(num & 1)
        num >>= 1
    return bits


def take_nth(ls, n):
    return list(map(lambda x: x[n], ls))


def read_traces(ntraces, nshares, file, file_bc, ignore_bc):
    traces = []
    with open(file, 'rb') as f:
        trace_num = 2*int.from_bytes(f.read(4), byteorder="big")
        samples_per_trace = int.from_bytes(f.read(4), byteorder="big")
    traces = np.fromfile(file, offset=8, dtype=np.float64)
    trace_num_2 = traces.shape[0]//samples_per_trace
    assert trace_num == trace_num_2
    assert trace_num == 2*ntraces, f"{trace_num} != {ntraces}"
    print(f"Read {traces.shape} floats.")
    print(f"{trace_num=}")
    traces = traces.reshape((trace_num, samples_per_trace))
    if ignore_bc:
        bcs = None
    else:
        bcs = np.fromfile(file_bc, offset=4, dtype=np.uint32,
                      count=SIMPLECOMPBITS*nshares*trace_num)
        bcs = bcs.reshape((trace_num, SIMPLECOMPBITS, nshares))
    return traces, bcs


def print_plt(x, y, name="plot", path="../plot_prints"):
    now = datetime.now()
    dt_string = now.strftime("%d%m%Y-%h%m%s")
    ls = list(zip(x, y))
    fname = f"{path}/{name}_{dt_string}.txt"
    print(f"Printing trace to {fname}..")
    try:
        with open(fname, 'w') as f:
            for i in range(len(ls)):
                if i % 5000 == 0:
                    f.write("\n")
                f.write(f"({ls[i][0]}, {ls[i][1]})")
    except OSError as e:
        print(f"Unable to write {fname}: {e}", file=sys.stderr)
        print("Are you sure the plot printing directory exists?")


def plot_traces(traces, titles, locs=None, sharey=False, sharex=False):
    if len(traces) == 1:
        fig, ax = plt.subplots(len(traces), sharey=sharey, sharex=sharex)
        plot_trace(traces[0], titles[0], ax, locs)
    else:
        fig, axs = plt.subplots(len(traces), sharey=sharey, sharex=sharex)
        for trace, title, ax in zip(traces, titles, axs):
            plot_trace(trace, title, ax, locs)
    plt.show()


def plot_trace(trace, title, ax, locs=None):
    if type(trace) is tuple:
        ax.plot(trace[0], trace[1])
        return
    else:
        ax.plot(trace)
    ax.set_title(title)
    print_plt(list(range(trace.shape[0])), trace, name=title.replace(" ", "_"))
    if locs is not None:
        for loc in locs:
            if type(loc) is not int:
                if loc.trace_loc is None:
                    continue
                loc = loc.trace_loc
            ax.plot([loc, loc], [max(trace), max(trace) +
                    0.1*(max(trace)-min(trace))], color="red")


def compute_auto_correlation(traces, max_len=1000, max_dist=5000, min_dist=100, start=300):
    cors_outer = []
    for loc0 in range(start, min(traces.shape[1], max_len)):
        cors = []
        for loc1 in range(loc0+min_dist, min(traces.shape[1], loc0+max_dist), min_dist):
            cor = np.cov(traces[:, loc0], traces[:, loc1])[0][1]
            cors.append((loc1, cor))
        cors = list(sorted(cors, key=lambda x: x[1], reverse=True))[:32]
        if len(cors) != 32:
            raise ValueError
        cors = list(sorted(cors, key=lambda x: x[0]))
        sum_cors = sum(c[1] for c in cors)
        cors_outer.append((loc0, cors, sum_cors))
    best_cors = list(sorted(cors_outer, key=lambda x: x[2], reverse=True))
    return best_cors


def separate_normals(samples):
    mean = np.mean(samples)
    indices_0 = np.where(samples >= mean)
    indices_1 = np.where(samples <= mean)
    mean_0 = np.mean(samples[indices_0])
    mean_1 = np.mean(samples[indices_1])
    sigma_0 = np.sqrt(np.var(samples[indices_0]))
    sigma_1 = np.sqrt(np.var(samples[indices_1]))
    return mean, (indices_0, mean_0, sigma_0), (indices_1, mean_1, sigma_1)


def find_poi_manual(traces, start=475, end=500, sel=False, shares=None, min_dist=50):
    from poi import Pois, Loc, PoisCollection
    print(f"Manual pois in range: {start}:{end}")
    mean = np.mean(traces, axis=0)
    if sel:
        assert shares is not None
        current_bit = 0
        current_share = 0
        pois = [[[None for _ in range(32)] for _ in range(shares)]]
        locs = []
        rem_locs = []
    loc = start
    while loc < min(traces.shape[1], end+1):
        samples = traces[:, loc]
        _, (idc0, _, _), (idc1, _, _) = separate_normals(samples)
        means0 = np.mean(traces[idc0], axis=0)
        means1 = np.mean(traces[idc1], axis=0)
        diff = means0 - means1
        print(f"Location in trace: {loc}")
        obs = traces[:, loc]
        hist0, bins0 = np.histogram(obs, density=True, bins=len(obs)//4)
        plot_traces([diff, (bins0[:-1], hist0), mean[max(0, loc-10*min_dist):min(mean.shape[0], loc+10*min_dist)], mean[0:end]], [f"Manual traces finding at {loc}", "Dist", "Mean local", "Mean with current loc"], locs=[loc])
        if sel:

            def print_cmd():
                inp = input(f"Select ({current_share}, {current_bit})? [y(es), (c)hange bit loc, (m)ove trace loc, (s)et skip, (r)emember, (p)rint]: ").strip()
                return inp

            inp = print_cmd()

            if inp == "m":
                inp = input("Set location in trace: ").strip()
                loc = int(inp)
                continue

            if inp == "c":
                inp = input("Set [current_share, current_bit]: ").strip()
                current_share = int(inp.split(',')[0].replace('(', "").strip())
                current_bit = int(inp.split(',')[1].replace(')', "").strip())
                inp = print_cmd()

            if inp == "s":
                inp = input("Set skip distance: ").strip()
                min_dist = int(inp)
                inp = print_cmd()

            if inp == "r":
                rem_locs.append(loc)
                inp = print_cmd()

            if inp == "p":
                print(f"POIs: {locs}")
                print(f"Remembered locations: {rem_locs}")
                inp = print_cmd()

            if inp == "y":
                print(f"Selected {loc} for {current_share=} {current_bit=}.")
                locs.append(loc)
                loc_cls = Loc(0, current_share, current_bit, loc)
                pois[0][current_share][current_bit] = Pois([loc_cls])
                current_bit += 1
                if current_bit == 32:
                    current_bit = 0
                    current_share += 1
                print(f"Selected until now: {locs}")
                print(f"Skipping next {min_dist} locations.")
                loc += min_dist
        loc += 1

    if sel:
        pois_col = PoisCollection(1, shares)
        pois_col.pois = pois
        return pois_col


def find_poi_separate_samples(traces):
    sep_trace_0 = [None for _ in range(traces.shape[1])]
    sep_trace_1 = [None for _ in range(traces.shape[1])]
    score_trace = [None for _ in range(traces.shape[1])]
    for loc in range(traces.shape[1]):
        samples = traces[:, loc]
        mean_0, mean_1, score = compute_separation_score(samples, loc)
        sep_trace_0[loc] = mean_0
        sep_trace_1[loc] = mean_1
        score_trace[loc] = score
    return np.array(sep_trace_0), np.array(sep_trace_1), np.array(score_trace)


def compute_separation_score(samples, loc):
    mean, (indices_0, mean_0, _), (indices_1, mean_1, _) = separate_normals(samples)
    score = (mean_1 - mean_0)**2
    l0 = indices_0[0].shape[0]
    l1 = indices_1[0].shape[0]
    rel = l0/(l0+l1)
    if rel < 0.45 or rel > 0.55:
        score = -0.01
    return mean_0, mean_1, score
