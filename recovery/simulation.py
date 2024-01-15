import numpy as np
from poi import Loc, Pois, PoisCollection
from util import bits, bits_2, share_value
from attack import attack


class Simulator:
    def __init__(self, sigma, num_bcs, num_shares, pois_per_bit):
        sigma_bound = 0.0001
        assert sigma >= sigma_bound, f"Sigma must be greater than {sigma_bound}"
        self.sigma = sigma
        self.num_bcs = num_bcs
        self.num_shares = num_shares
        self.pois_per_bit = pois_per_bit
        self.total_points = num_bcs*num_shares*pois_per_bit*32

        self.bcs = []
        self.traces = []
        self.pois = None

    def find_pois(self):
        assert self.pois is None
        self.pois = PoisCollection(self.num_bcs, self.num_shares)
        trace_loc = 0
        for bc_idx in range(self.num_bcs):
            for share_idx in range(self.num_shares):
                for bit_idx in range(32):
                    locs = []
                    for poi_idx in range(self.pois_per_bit):
                        loc = Loc(bc_idx, share_idx, bit_idx, trace_loc=trace_loc)
                        locs.append(loc)
                        trace_loc += 1
                    poi = Pois(locs)
                    self.pois.set_pois(bc_idx, share_idx, bit_idx, poi)
        assert self.pois.len() == self.num_bcs*self.num_shares*32
        assert self.pois.len() == self.total_points//self.pois_per_bit

    def record_trace(self, dec_fail=False, append=True):
        trace = []
        bcs = []
        for bc_idx in range(self.num_bcs):
            if dec_fail:
                shares = [np.random.randint(0, 1 << 32) for _ in range(self.num_shares)]
            else:
                shares = share_value(0, self.num_shares)
            bcs.append(shares)
            r = np.random.randint(0, 1 << 64, dtype=np.uint64)
            hw = sum(bits_2(r))
            for share_idx, share in enumerate(shares):
                current_share_bits = bits(share)
                for bit_idx in range(32):
                    hw_cur = 0 if current_share_bits[bit_idx] == 0 else hw
                    for poi_idx in range(self.pois_per_bit):
                        obs = np.random.normal(hw_cur, self.sigma)
                        trace.append(obs)
        assert len(trace) == self.total_points
        if append:
            self.traces.append(trace)
            self.bcs.append(bcs)
        return trace, bcs

    def record_traces(self, num_traces, verbose=True):
        for i in range(num_traces):
            self.record_trace(i & 1 == 1)
            if verbose:
                print(f"{i}/{num_traces}" + " "*20, end='\r')

    def finish_recording_phase(self):
        self.traces = np.array(self.traces)
        self.bcs = np.array(self.bcs)

    def execute_attacks(self, number_of_traces=None, profile_traces=None, template=True, vertical=True, horizontal=True):
        res = {"sigma": self.sigma, "traces": len(self.traces), "nshares": self.num_shares, "nbcs": self.num_bcs, "npois": self.pois_per_bit}
        if template:
            res["results_template"] = self.execute_attack("template", number_of_traces=number_of_traces, profile_traces=profile_traces)
        if vertical:
            res["results_vertical"] = self.execute_attack("auto_vertical", number_of_traces=number_of_traces)
        if horizontal:
            res["results_horizontal"] = self.execute_attack("auto_horizontal", number_of_traces=number_of_traces)
        return res

    def execute_attack(self, attack_name, profile_traces=None, reset=True, number_of_traces=None):
        if attack_name == "template":
            print("######## SIMULATED TEMPLATE ATTACK ########")
            print(f"Sampling {2*profile_traces} template traces..")
            template_traces = []
            template_bcs = []
            for _ in range(profile_traces):
                trace0, bc0 = self.record_trace(False, append=False)
                trace1, bc1 = self.record_trace(True, append=False)
                template_traces.append(trace0)
                template_traces.append(trace1)
                template_bcs.append(bc0)
                template_bcs.append(bc1)
            template_traces = np.array(template_traces)
            template_bcs = np.array(template_bcs)

            def tmpl_func():
                return self.pois.compute_template_from_bc(template_traces, template_bcs)
            total, correct, ratio_correct, classified, ratio_classified = attack(self.traces, self.pois, self.num_shares, f"Simulated templated attack {self.num_shares}-{len(self.traces)}-{self.sigma}", template_function=tmpl_func, number_of_traces=number_of_traces)
        elif attack_name == "auto_vertical":
            def tmpl_func():
                return self.pois.compute_vertical_auto_template(self.traces)
            total, correct, ratio_correct, classified, ratio_classified = attack(self.traces, self.pois, self.num_shares, f"Simulated vertical auto template attack {self.num_shares}-{len(self.traces)}-{self.sigma}", template_function=tmpl_func, number_of_traces=number_of_traces)
        elif attack_name == "auto_horizontal":
            total, correct, ratio_correct, classified, ratio_classified = attack(self.traces, self.pois, self.num_shares, f"Simulated horizontal auto template attack {self.num_shares}-{len(self.traces)}-{self.sigma}", build_single_trace_template=PoisCollection.compute_horizontal_auto_template, number_of_traces=number_of_traces, template_function=None)
        else:
            raise ValueError("Unknown attack")
        if reset:
            self.pois.reset_template()
        return total, correct, ratio_correct, classified, ratio_classified
