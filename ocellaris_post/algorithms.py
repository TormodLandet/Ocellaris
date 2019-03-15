# Copyright (C) 2018-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0

def get_stairs(all_results, report_name, min_length=50):
    """
    Given a list of Result objects and the name of a control variable
    time series, return two lists

    - all_vals: all discovered flat plateau values in the report_name time
      series. The plateau valuse for all reults are merged into one ordered
      list with unique values only

    - all_lims: same length as all_vals; each item is a list of the same
      length as all_results and each list inside this inner list is a tuple of
      start and end indices of the flat part; or (0, 0, 0) if there was no such
      flat part in that result object

    The indices returned are (istart_trans, istart, iend):

    - istart_trans: where the unsteady transition starts ("step up")
    - istart: where the steady part of this "stair tread" starts
    - iend: where the steady part of this "stair tread" ends
    """
    # Find start and end indices of flat parts of the control value time series
    all_res = []
    all_vals = set()

    for results in all_results:
        this_res = []
        all_res.append(this_res)

        ts = results.reports.get(report_name, [])
        N = len(ts)

        start = True
        for i in range(N):
            if start and i < N - 2 and ts[i] == ts[i + 1] == ts[i + 2]:
                start = False
                this_res.append([ts[i], i, None])
                all_vals.add(ts[i])
                continue
            elif (
                not start
                and N - 1 > i > 1
                and ts[i - 2] == ts[i - 1] == ts[i] != ts[i + 1]
            ):
                start = True
                this_res[-1][-1] = i

        if not start:
            this_res[-1][-1] = N - 1

    # For each flat control value, make lists of start and stop indices for each result object
    all_vals = sorted(all_vals)
    all_lims = []
    for val in all_vals:
        all_lims.append([])

        for lims in all_res:
            prevend = 0
            for v, istart, iend in lims:
                if v == val and iend - istart > min_length:
                    all_lims[-1].append((prevend, istart, iend))
                    break
                prevend = iend
            else:
                all_lims[-1].append((0, 0, 0))

    # Remove values that ended up with no data due to min_length
    all_vals2, all_lims2 = [], []
    for val, lims in zip(all_vals, all_lims):
        all_zero = True
        for i0, i1, i2 in lims:
            if not (i0 == i1 == i2 == 0):
                all_zero = False
        if not all_zero:
            all_vals2.append(val)
            all_lims2.append(lims)

    return all_vals2, all_lims2
