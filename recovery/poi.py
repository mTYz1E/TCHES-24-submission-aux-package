import copy
import numpy as np
import scipy
from util import bits, take_nth, separate_normals


class Loc:
    def __init__(self, bc_index, share_index, bit_index, trace_loc=None):
        self.bc_index = bc_index
        self.share_index = share_index
        self.bit_index = bit_index
        self.trace_loc = trace_loc

    def __str__(self):
        return f"{self.bc_index=}, {self.share_index=}, {self.bit_index=}, {self.trace_loc=}"

    def get_by_bc(self, traces, bcs):
        bit = self.bit_index
        bc = np.array([bits(bc_e)[bit] for bc_e in bcs[:, self.bc_index, self.share_index]])
        indices_0 = np.argwhere(bc == 0)[:, 0]
        indices_1 = np.argwhere(bc == 1)[:, 0]
        assert indices_0.shape[0] + indices_1.shape[0] == traces.shape[0]
        return indices_0, indices_1

    def compute_mean_and_var(self, traces, bcs):
        indices_0, indices_1 = self.get_by_bc(traces, bcs)
        assert indices_0.shape[0] + indices_1.shape[0] == traces.shape[0]
        means_0 = np.mean(traces[indices_0], axis=0)
        means_1 = np.mean(traces[indices_1], axis=0)
        if self.trace_loc is not None:
            mean_0, mean_1 = means_0[self.trace_loc], means_1[self.trace_loc]
            var_0 = np.var(traces[indices_0][:, self.trace_loc])
            var_1 = np.var(traces[indices_1][:, self.trace_loc])
            return means_0, means_1, mean_0, mean_1, var_0, var_1
        return means_0, means_1, None, None, None, None


class Pois:
    def __init__(self, locs):
        assert len(locs) >= 1
        for loc in locs:
            assert loc.bc_index == locs[0].bc_index
            assert loc.share_index == locs[0].share_index
            assert loc.bit_index == locs[0].bit_index
        self.locs = locs
        self.trace_locs = np.array([lo.trace_loc for lo in locs])
        self.template = None

    def compute_vertical_auto_template(self, traces):
        # Separate the distributions and classify traces according to it
        indices_0 = []
        indices_1 = []
        separated_loc_dists = []
        for loc in self.trace_locs:
            samples = traces[:, loc]
            mean, (_, mean_0, sigma_0), (_, mean_1, sigma_1) = separate_normals(samples)
            separated_loc_dists.append((mean, (mean_0, sigma_0), (mean_1, sigma_1)))

        for i, trace_i in enumerate(traces):
            ev0 = 1
            ev1 = 1
            diffs = []
            for loc, (mean, (mean_0, sigma_0), (mean_1, sigma_1)) in zip(self.trace_locs, separated_loc_dists):
                obs = trace_i[loc]
                # Does not work as 1-bits are not normally distributed
                # v0 = scipy.stats.norm(mean_0, sigma_0).pdf(obs)
                # v1 = scipy.stats.norm(mean_1, sigma_1).pdf(obs)
                v0 = 2 if obs > mean else 0
                v1 = 2 if obs < mean else 0
                # v0 = 2 if obs >= mean else 0
                # v1 = 2 if obs <= mean else 0
                ev0 *= v0
                ev1 *= v1
                diffs.append((abs(obs-mean_0), abs(obs-mean_1)))
            if ev0 >= 2:
                indices_0.append(i)
            elif ev1 >= 2:
                indices_1.append(i)
            else:
                if diffs[0] < diffs[1]:
                    indices_0.append(i)
                else:
                    indices_1.append(i)

        self.compute_template(traces, indices_0, indices_1)

    def compute_template_from_bc(self, traces, bcs):
        assert self.template is None
        if self.get_num_pois() == 1:
            _, _, mean_0, mean_1, var_0, var_1 = self.locs[0].compute_mean_and_var(traces, bcs)
            self.template = (([mean_0], [var_0]), ([mean_1], [var_1]))
        else:
            indices_0, indices_1 = self.locs[0].get_by_bc(traces, bcs)
            self.compute_template(traces, indices_0, indices_1)

    def compute_template(self, traces, indices_0, indices_1):
        assert self.template is None
        means_0 = []
        means_1 = []
        cov_0 = [[None for _ in range(self.get_num_pois())] for _ in range(self.get_num_pois())]
        cov_1 = [[None for _ in range(self.get_num_pois())] for _ in range(self.get_num_pois())]
        for i, loc_a in enumerate(self.locs):
            mean_a_0 = np.mean(traces[indices_0, loc_a.trace_loc])
            mean_a_1 = np.mean(traces[indices_1, loc_a.trace_loc])
            means_0.append(mean_a_0)
            means_1.append(mean_a_1)
            for j, loc_b in enumerate(self.locs):
                cov_matrix_0 = np.cov(traces[indices_0, loc_a.trace_loc], traces[indices_0, loc_b.trace_loc])
                cov_matrix_1 = np.cov(traces[indices_1, loc_a.trace_loc], traces[indices_1, loc_b.trace_loc])
                cov_0[i][j] = cov_matrix_0[0][1]
                cov_1[i][j] = cov_matrix_1[0][1]
        self.template = ((means_0, cov_0), (means_1, cov_1))

    def apply_template(self, trace):
        assert self.template is not None
        obs = trace[self.trace_locs]
        if self.get_num_pois() == 1:
            mean_0, sig_0 = self.template[0][0], self.template[0][1]
            mean_1, sig_1 = self.template[1][0], self.template[1][1]
            v0 = scipy.stats.norm(mean_0[0], np.sqrt(sig_0[0])).pdf(obs[0])
            v1 = scipy.stats.norm(mean_1[0], np.sqrt(sig_1[0])).pdf(obs[0])
        else:
            mean_0, cov_0 = self.template[0][0], self.template[0][1]
            mean_1, cov_1 = self.template[1][0], self.template[1][1]
            v0 = scipy.stats.multivariate_normal(mean_0, cov_0).pdf(obs)
            v1 = scipy.stats.multivariate_normal(mean_1, cov_1).pdf(obs)
        bit = 1 if v1 > v0 else 0
        return bit, v0, v1

    def __str__(self):
        return f"POI at {self.locs}"

    def get_num_pois(self):
        return len(self.locs)

    def reset_template(self):
        self.template = None

    def set_template(self, template):
        self.template = template

    @classmethod
    def find_pois(cls, traces, bcs, loc, num_pois=1, min_distance=5):
        means0, means1, _, _, _, _ = loc.compute_mean_and_var(traces, bcs)
        diff = means0 - means1
        sort = np.argsort(np.abs(diff), axis=0)
        found = 0
        locs = []
        for i in range(1, min_distance+1):
            x = sort[-i]
            for p in locs:
                if abs(x - p.trace_loc) < min_distance:
                    continue
            loc_cur = copy.deepcopy(loc)
            loc_cur.trace_loc = x
            locs.append(loc_cur)
            found += 1
            if found >= num_pois:
                break
        if found < num_pois:
            print("Did not find sufficiently many POIs")
        return cls(locs)


class PoisCollection:
    def __init__(self, num_bcs, shares):
        self.pois = [[[None for _ in range(32)] for _ in range(shares)] for _ in range(num_bcs)]

    def len(self):
        i = 0

        def count(_):
            nonlocal i
            i += 1
        self.map(count)
        return i

    def get_pois_list(self):
        ls = [pbit for pbc in self.pois for pshares in pbc for pbit in pshares]
        return ls

    def get_num_bcs(self):
        return len(self.pois)

    def get_num_shares(self):
        return len(self.pois[0])

    def get_num_pois_per_bit(self):
        return self.pois[0][0][0].get_num_pois()

    def get_poi(self, idx_bc, share, bit):
        return self.pois[idx_bc][share][bit]

    def set_pois(self, idx_bc, share, bit, poi):
        self.pois[idx_bc][share][bit] = poi

    def reset_template(self):
        self.map(Pois.reset_template)

    def compute_template_from_bc(self, traces, bcs):
        self.map(Pois.compute_template_from_bc, traces, bcs)

    def compute_vertical_auto_template(self, traces):
        self.map(Pois.compute_vertical_auto_template, traces)

    def compute_auto_template_simple(self, traces):
        self.map(Pois.compute_auto_template_simple, traces)

    def compute_horizontal_auto_template(self, trace):
        for bc_idx in range(self.get_num_bcs()):
            self.compute_horizontal_auto_template_bc(trace, bc_idx)

    def compute_horizontal_auto_template_bc(self, trace, bc_idx):
        num_pois = self.get_num_pois_per_bit()
        pois = [poi_bit for pois_share in self.pois[bc_idx] for poi_bit in pois_share]
        locs = [[p.locs[poi_idx].trace_loc for p in pois] for poi_idx in range(num_pois)]
        obs = np.array([[trace[loc] for loc in locs[poi_idx]] for poi_idx in range(num_pois)])

        separated_loc_dists = []

        for poi_idx in range(num_pois):
            samples = obs[poi_idx]
            mean, (_, mean_0, sigma_0), (_, mean_1, sigma_1) = separate_normals(samples)
            separated_loc_dists.append((mean, (mean_0, sigma_0), (mean_1, sigma_1)))

        indices_0 = []
        indices_1 = []

        for i in range(len(pois)):
            ev0 = 1
            ev1 = 1
            diffs = []
            for poi_idx, (mean, (mean_0, sigma_0), (mean_1, sigma_1)) in enumerate(separated_loc_dists):
                ob = obs[poi_idx][i]
                # Does not work as 1-bits are not normally distributed
                # v0 = scipy.stats.norm(mean_0, sigma_0).pdf(obs)
                # v1 = scipy.stats.norm(mean_1, sigma_1).pdf(obs)
                v0 = 2 if ob > mean else 0
                v1 = 2 if ob < mean else 0
                ev0 *= v0
                ev1 *= v1
                diffs.append((abs(ob-mean_0), abs(ob-mean_1)))
            if ev0 >= 2:
                indices_0.append(i)
            elif ev1 >= 2:
                indices_1.append(i)
            else:
                if sum(take_nth(diffs, 0)) < sum(take_nth(diffs, 1)):
                    indices_0.append(i)
                else:
                    indices_1.append(i)

        obs_0 = obs[:, indices_0]
        obs_1 = obs[:, indices_1]
        means_0 = [None for _ in range(num_pois)]
        means_1 = [None for _ in range(num_pois)]
        covs_0 = [[None for _ in range(num_pois)] for _ in range(num_pois)]
        covs_1 = [[None for _ in range(num_pois)] for _ in range(num_pois)]
        for i in range(num_pois):
            mean_0 = np.mean(obs_0[i])
            mean_1 = np.mean(obs_1[i])
            means_0[i] = mean_0
            means_1[i] = mean_1
            for j in range(num_pois):
                obs_i_0 = obs_0[i]
                obs_i_1 = obs_1[i]
                obs_j_0 = obs_0[j]
                obs_j_1 = obs_1[j]
                cov_0 = np.cov(obs_i_0, obs_j_0)[0][1]
                cov_1 = np.cov(obs_i_1, obs_j_1)[0][1]
                covs_0[i][j] = cov_0
                covs_1[i][j] = cov_1
        for poi in pois:
            poi.set_template(((means_0, np.array(covs_0)), (means_1, np.array(covs_1))))

    def apply_template(self, trace):
        return self.map(Pois.apply_template, trace)

    def map(self, fn, *args, **kargs):
        results = []
        for bc_pois in self.pois:
            results_bc = []
            for share_pois in bc_pois:
                results_share = []
                for bit_pois in share_pois:
                    res = fn(bit_pois, *args, **kargs)
                    results_share.append(res)
                results_bc.append(results_share)
            results.append(results_bc)
        return results

    @classmethod
    def find_all_pois(cls, traces, bcs, num_bcs=1, num_pois_per_bit=1, min_distance=5):
        pois = cls(num_bcs, bcs.shape[2])
        for idx_bc in range(num_bcs):
            for share in range(bcs.shape[2]):
                for bit_idx in range(32):
                    loc = Loc(idx_bc, share, bit_idx)
                    pois_loc = Pois.find_pois(traces, bcs, loc, num_pois_per_bit, min_distance)
                    pois.set_pois(idx_bc, share, bit_idx, pois_loc)
        return pois
