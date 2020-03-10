import aposteriori_error
import numpy
from pathlib import Path, PurePath
import glob
import pprint
import matplotlib.pyplot as plt
import json


def read_time():
    with open("../training/validation/time_trajectory_39/time.dat") as fi:
        for line in fi.readlines():
            if "User time" in line:
                utime = float(line.split(":")[-1])
                continue
            if "System time" in line:
                stime = float(line.split(":")[-1])
                continue
        time_hf = utime + stime

    case_paths = glob.glob("../multiscale_1ip_speedup/case_*/time.dat")
    times = {}
    for case_path in case_paths:
        with open(case_path) as fi:
            for line in fi.readlines():
                if "User time" in line:
                    utime = float(line.split(":")[-1])
                    continue
                if "System time" in line:
                    stime = float(line.split(":")[-1])
                    continue
        case = case_path.split("/")[2]
        times[case] = utime + stime
    return times, time_hf


def read_data():
    stress_ref = numpy.loadtxt(
        "../training/validation/trajectory_39/homogenized_stress.dat"
    )
    case_paths = glob.glob("../multiscale_1ip/case_*/homogenized_stress.dat")
    errors = {}
    for case_path in case_paths:
        stress = numpy.loadtxt(case_path)
        case = case_path.split("/")[2]
        errors[case] = max(
            aposteriori_error.compute_mean_normalized_error(stress, stress_ref)[6:]
        )
    return errors


def make_matrix(errors):
    m = []
    p = []
    for k, v in errors.items():
        mi = int(k.split("_")[2].split("m")[0])
        m.append(mi)
        pi = int(k.split("_")[3].split("ip")[0])
        p.append(pi)
    ms = list(sorted(set(m)))
    ps = list(sorted(set(p)))
    a = {}
    for m in ms:
        a[m] = []
        for p in ps:
            k = "case_39t_{}m_{}ip".format(m, p)
            a[m].append(errors[k])
    return a, ps


##################################
# main
##################################
if __name__ == "__main__":
    fp = open("config.json", "r")
    params = json.load(fp)
    errors = read_data()
    times, time_hf = read_time()
    print("Time HF: {}s".format(time_hf))
    first = True
    for c in sorted(errors, key=errors.get):
        if first:
            er = errors[c]
            tr = times[c]
            first = False
            print(
                "{:<25} {}      {}      {}  {}  {}  {}".format(
                    "Case", "error", "time", "s/err", "s=tr/t", "err=e/er", "speedup"
                )
            )
            print("----------------------------------------------------------------")
        e = errors[c]
        t = times[c]
        s = tr / t
        err = e / er
        print(
                "{:<25} {:>1.3e} {:>6.2f}s {:>8.3f} {:>8.3f} {:>8.3f} {:>8.3f}".format(
                c, e, t, s / err, s, err, time_hf/t
            )
        )
    # errors_series = make_matrix(errors)
    a, p = make_matrix(errors)
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=False)
    for k, v in a.items():
        ax1.plot(p, v, label=k, marker="")
        ax2.plot(p, v, label=k, marker="o")
    # configure plot
    ax1.legend()
    #ax1.set_title(params["ax1_title"])
    #ax2.set_ylabel("error")
    ax2.set_ylim(0, params["ax2_ylimit"])
    # ax1.spines['bottom'].set_visible(False)
    # ax2.spines['top'].set_visible(False)
    # ax1.xaxis.tick_top()
    # ax1.tick_params(labeltop=False)  # don't put tick labels at the top
    # ax2.xaxis.tick_bottom()

    plt.show()
