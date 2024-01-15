import numpy as np
import scipy
import matplotlib.pyplot as plt
from util import bits, bits_2, print_plt, plot_traces, compute_auto_correlation, take_nth, find_poi_separate_samples
from poi import Loc


def plot_distribution_bc_horizontal(trace, pois, bcs=None, idx_bc=0, idx_poi=0, inverse_bin_fac=2):
    locs = []
    locs_1 = []
    shares = pois.get_num_shares()
    for share in range(shares):
        for bit in range(32):
            loc = pois.get_poi(idx_bc, share, bit).locs[idx_poi].trace_loc
            if bcs is not None:
                b = bits(bcs[idx_bc][share])[bit]
                if b == 0:
                    locs.append(loc)
                else:
                    locs_1.append(loc)
            else:
                locs.append(loc)
    if bcs is None:
        obs = trace[locs]
        plt.hist(obs, bins=len(obs))
    else:
        obs_0 = trace[locs]
        obs_1 = trace[locs_1]
        # plt.hist(obs_0, bins=len(obs_0))
        # plt.hist(obs_1, bins=len(obs_1))
        hist0, bins0 = np.histogram(obs_0, density=False, bins=len(obs_0)//inverse_bin_fac)
        hist1, bins1 = np.histogram(obs_1, density=False, bins=len(obs_1)//inverse_bin_fac)
        plt.plot(bins0[:-1], hist0)
        plt.plot(bins1[:-1], hist1)
        print_plt(bins0, hist0, name="horizontal_0")
        print()
        print_plt(bins1, hist1, name="horizontal_1")
    plt.title(f"Horizontal distribution with {len(locs) + len(locs_1)} samples")
    plt.show()


def plot_distribution_bc_vertical(traces, pois, bcs=None, trace_idx=0, idx_bc=0, share=0, bit=0, idx_poi=0, inverse_bin_fac=4, plot_loc=None):
    if bcs is None:
        if plot_loc is None:
            loc = pois.get_poi(idx_bc, share, bit).locs[idx_poi].trace_loc
        else:
            loc = plot_loc
        obs = traces[:, loc]
        plt.hist(obs, bins=len(obs)//inverse_bin_fac)
    else:
        loc = pois.get_poi(idx_bc, share, bit).locs[idx_poi]
        mean = np.mean(traces[:, loc.trace_loc])
        idx0, idx1 = loc.get_by_bc(traces, bcs)
        obs_0 = traces[idx0, loc.trace_loc]
        obs_1 = traces[idx1, loc.trace_loc]
        incorrect_0 = obs_0[obs_0 < mean]
        incorrect_1 = obs_1[obs_1 > mean]
        mean_0 = np.mean(obs_0)
        mean_1 = np.mean(obs_1)
        sigma_0 = np.sqrt(np.var(obs_0))
        sigma_1 = np.sqrt(np.var(obs_1))
        print(f"Mean of all traces at the POI: {mean}")
        print(f"Number of zero bits below mean: {incorrect_0.shape[0]}")
        print(f"Number of one bits above mean: {incorrect_1.shape[0]}")
        print(f"{mean_0=}, {mean_1=}, {sigma_0=}, {sigma_1=}")
        # plt.hist(obs_0, bins=len(obs_0)//2)
        # plt.hist(obs_1, bins=len(obs_1)//2)
        hist0, bins0 = np.histogram(obs_0, density=True, bins=len(obs_0)//inverse_bin_fac)
        hist1, bins1 = np.histogram(obs_1, density=True, bins=len(obs_1)//inverse_bin_fac)
        plt.plot(bins0[:-1], hist0)
        plt.plot(bins1[:-1], hist1)
        print_plt(bins0, hist0, name="vertical_0.txt")
        print_plt(bins1, hist1, name="vertical_1.txt")
        # plt.plot([mean, mean], [0, 1])
    plt.title(f"Vertical distribution with {len(traces)} samples")
    plt.show()


def plot_modeled_vertical_distributions(num, sigma):
    obss = []
    rs = []
    hws = []
    obss_0 = []
    for i in range(num):
        r = np.random.randint(0, 1 << 64, dtype=np.uint64)
        rs.append(r)
        hw = sum(bits_2(r))
        hws.append(hw)
        obs = np.random.normal(hw, sigma)
        obss.append(obs)

        obss_0.append(np.random.normal(0, sigma))
    sigma = np.sqrt(np.var(obss))
    mu = np.mean(obss)
    print(f"{mu=}, {sigma=}")
    hist, bins = np.histogram(obss, bins=list(range(-32, 64)), density=True)
    hist0, bins0 = np.histogram(obss_0, bins=list(range(-32, 64)), density=True)
    x = np.linspace(0, 64, 100)
    y = scipy.stats.norm.pdf(x, mu, sigma)
    plt.plot(x, y)
    print_plt(bins, hist, name="sim_vertical.txt")
    print_plt(x, y)
    plt.plot(bins[:-1], hist)
    print_plt(bins0, hist0)
    plt.plot(bins0[:-1], hist0)
    plt.show()


def plot_diff_means(traces, bcs, bc_idx=0, share=0, bit=0, pois=None):
    if bcs is not None:
        poi = None
        if pois is not None:
            poi = pois.get_poi(bc_idx, share, bit).locs
        loc = Loc(bc_idx, share, bit)
        m0, m1, _, _, _, _ = loc.compute_mean_and_var(traces, bcs)
        plot_traces([np.mean(traces, axis=0), m0-m1], ["Mean", f"Difference of Means {bc_idx=} {share=} {bit=}"], locs=poi)
    else:
        poi = None
        if pois is not None:
            poi = pois.get_poi(bc_idx, share, bit).locs
        plot_traces([np.mean(traces, axis=0)], ["Mean"], locs=poi)


def plot_auto_correlation(traces, bcs, pois):
    cors = compute_auto_correlation(traces)
    loc0, others0, _ = cors[0]
    locs_raw = [loc0] + take_nth(others0, 0)
    locs = list(map(lambda x: Loc(0, 0, 0, x), locs_raw))
    loc = Loc(0, 0, 0)
    m0, m1, _, _, _, _ = loc.compute_mean_and_var(traces, bcs)
    print(list(map(lambda x: x.trace_locs[0], pois)))
    print(locs_raw)
    plot_traces([m0-m1], ["Difference of Means"], locs=locs)


def plot_separated_pois(traces, bcs=None, pois=None):
    sep0, sep1, score = find_poi_separate_samples(traces)
    locs_greater_0 = np.where(score > 0)
    print(locs_greater_0)
    print(list(map(lambda x: x.trace_locs[0], pois)))
    loc = Loc(0, 0, 0)
    m0, m1, _, _, _, _ = loc.compute_mean_and_var(traces, bcs)
    plot_traces([m0-m1, sep0, sep1, score], ["Difference of Means", sep0, sep1, score])


def plot_t_test_fail_nfail(traces, pois=None):
    idc_0 = list(range(0, traces.shape[0], 2))
    idc_1 = list(range(1, traces.shape[0], 2))
    plot_t_test(traces, idc_0, idc_1, pois)


def plot_t_test(traces, idc_0, idc_1, pois=None):
    t_0 = traces[idc_0, :]
    t_1 = traces[idc_1, :]
    mean_sucs = np.mean(t_0, axis=0)
    mean_fail = np.mean(t_1, axis=0)
    mean_diff = mean_sucs-mean_fail
    assert t_0.shape[0] + t_1.shape[0] == traces.shape[0]
    corrected_var_sucs = np.var(t_0, axis=0)/t_0.shape[0]
    corrected_var_fail = np.var(t_1, axis=0)/t_1.shape[0]
    t_stat = mean_diff/np.sqrt(corrected_var_fail+corrected_var_sucs)
    leakage_points = np.where(t_stat >= 4.5)[0]
    print(f"{leakage_points=}")
    if pois is not None:
        for poi_b in pois.get_pois_list():
            for poi in poi_b.locs:
                if poi.trace_loc in leakage_points:
                    print(f"WARNING: {poi} is at a leaking point!")

    plot_traces([np.mean(traces, axis=0), t_stat], ["mean", "t-statistic"], sharey=False, sharex=True)
