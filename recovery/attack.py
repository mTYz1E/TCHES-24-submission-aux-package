import numpy as np
from util import take_nth, bits


def eval_attack(res_is_success, name):
    correct = 0
    total = 0
    classified = 0
    unclassified = 0

    for i, suc in enumerate(res_is_success):
        total += 1
        if suc is None:
            unclassified += 1
            continue
        classified += 1
        if (i % 2 == 0 and suc) or (i % 2 == 1 and not suc):
            correct += 1
    assert classified + unclassified == total
    assert total == len(res_is_success)

    ratio_correct = correct/classified if classified != 0 else 0
    unclassified = total - classified
    ratio_classified = classified/total
    print(" "*40)
    print(f"{name}: {total=}, {classified=}, {unclassified=} {correct=}, {ratio_correct=}, {ratio_classified=}")
    return total, correct, ratio_correct, classified, ratio_classified


def attack(traces, pois, nshares, name="", build_single_trace_template=None, upper_bound=0.55, verbose=True, number_of_traces=None, test_mode=True, template_function=None):
    print(f"######## {name.upper()} ########")
    if template_function is not None:
        print("Building templates..")
        template_function()
    print("Applying templates to every trace..")
    res_is_success = []
    if number_of_traces is None:
        number_of_traces = traces.shape[0]
    else:
        number_of_traces = min(traces.shape[0], number_of_traces)
    for idx_tr in range(traces.shape[0]):
        if number_of_traces is not None and idx_tr >= number_of_traces:
            break
        if build_single_trace_template is not None:
            build_single_trace_template(pois, traces[idx_tr])
        rec = pois.apply_template(traces[idx_tr])
        if build_single_trace_template is not None:
            pois.reset_template()
        suc = classify_from_recovered(rec, idx_tr, upper_bound, nshares)
        res_is_success.append(suc)
        print(f"{idx_tr}/{number_of_traces}" + " "*20, end='\r')

    if test_mode:
        retv = eval_attack(res_is_success, name)
    else:
        print(res_is_success)
        retv = res_is_success

    print("Resetting template..")
    pois.reset_template()
    print("##############################################")
    print()

    return retv


def classify_from_recovered(rec, idx_tr, upper_bound, nshares):
    total = 0
    zeros = 0
    for idx_bc in range(len(rec)):
        final = np.array([0 for _ in range(32)])
        for idx_share in range(nshares):
            share = take_nth(rec[idx_bc][idx_share], 0)
            final = final + np.array(share)
        final = [b % 2 for b in final]
        zeros += len(final) - sum(final)
        total += len(final)
    ratio_correct = zeros/total
    if ratio_correct > 0.8:
        return True
        if idx_tr % 2 != 0:
            print(f"False positive: {idx_tr=}, {ratio_correct=}")
    elif ratio_correct < upper_bound:
        return False
        if idx_tr % 2 != 1:
            print(f"False negative: {idx_tr=}, {ratio_correct=}")
    else:
        return None


def attack_one_trace(traces, bcs, tr_idx, pois, name="", template_function=None):
    print(f"######## {name.upper()} ########")
    if template_function is not None:
        print("Building templates..")
        template_function()
    print("Applying template to trace..")
    rec = pois.apply_template(traces[tr_idx])
    compare_recovered(rec, bcs, tr_idx, name)
    print()
    print("Resetting template..")
    pois.reset_template()
    print("##############################################")
    print()


def compare_recovered(bcs_rec, bcs, tr_idx, name=""):
    correct = 0
    incorrect = 0
    unclassified = 0
    total = 0
    for idx_bc, bc_rec in enumerate(bcs_rec):
        for idx_share, bc_share in enumerate(bc_rec):
            bc_exp_cur = bits(bcs[tr_idx, idx_bc, idx_share])
            for bit_idx, (br, be) in enumerate(zip(bc_share, bc_exp_cur)):
                if br[0] is None:
                    unclassified += 1
                elif br[0] == be:
                    correct += 1
                else:
                    incorrect += 1
                    print(f"{name}: Incorrect bit at {idx_bc=}, {idx_share=}, {bit_idx=}, actual={be}, ev0={br[1]}, ev1={br[2]}")
                total += 1
    assert total == correct + incorrect + unclassified
    classified = total - unclassified
    ratio_correct = correct/total
    print(f"{name}: {total=}, {classified=}, {unclassified=}, {correct=}, {incorrect=}, {ratio_correct=}")
