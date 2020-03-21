import pprint
import aposteriori_error
import numpy
from pathlib import Path, PurePath
import glob
import pprint
import matplotlib.pyplot as plt
import json
import matplotlib.ticker as ticker


def read_time():
    with open("../training/validation/time_trajectory_99/time.dat") as fi:
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
        if "ROM" in case_path:
            continue
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
        "../training/validation/trajectory_99/homogenized_stress.dat"
    )
    case_paths = glob.glob("../multiscale_1ip/case_*/homogenized_stress.dat")
    errors = {}
    for case_path in case_paths:
        if "ROM" in case_path:
            continue
        stress = numpy.loadtxt(case_path)
        case = case_path.split("/")[2]
        errors[case] = max(
            aposteriori_error.compute_mean_normalized_error(stress, stress_ref)[6:]
        )
    return errors


def read_rom():
    stress_ref = numpy.loadtxt(
        "../training/validation/trajectory_99/homogenized_stress.dat"
    )
    case_paths = glob.glob("../multiscale_1ip/ROM_cases/case_*/homogenized_stress.dat")
    errors = {}
    for case_path in case_paths:
        stress = numpy.loadtxt(case_path)
        case = case_path.split("/")[3]
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
            k = "case_99t_{}m_{}ip".format(m, p)
            a[m].append(errors[k])
    return a, ps


##################################
# main
##################################
if __name__ == "__main__":
    fp = open("config.json", "r")
    params = json.load(fp)
    fo = open(params["dat_fname"], "w")
    errors = read_data()
    errors_rom = read_rom()
    times, time_hf = read_time()
    line = "Time HF: {}s".format(time_hf)
    print(line)
    fo.write(line + "\n")
    first = True
    for c in sorted(errors, key=errors.get):
        if first:
            er = errors[c]
            tr = times[c]
            first = False
            line = "{:<25} {}      {}      {}  {}  {}  {}".format( "case", "error", "time", "s/err", "s=tr/t", "err=e/er", "speedup")
            print(line)
            fo.write(line + "\n")
            line = "----------------------------------------------------------------"
            print(line)
            fo.write(line + "\n")
        e = errors[c]
        t = times[c]
        s = tr / t
        err = e / er
        line = "{:<25} {:>1.3e} {:>6.2f}s {:>8.3f} {:>8.3f} {:>8.3f} {:>8.3f}".format( c, e, t, s / err, s, err, time_hf/t)
        print(line)
        fo.write(line + "\n")
    # errors_series = make_matrix(errors)

    a, p = make_matrix(errors)
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=False)
    ax3 = ax2.twinx()
    for k, v in a.items():
        ax1.plot(p, v, label=k, marker="")
        ax2.plot(p, v, label=k, marker="o")
    color = 0
    error_rom = []
    for k, v in errors_rom.items():
        ax3.plot(p, [v] * len(p), label=k, marker="", linestyle='dotted', color="C{}".format(color))
        error_rom.append(v)
        color += 1

    # configure plot
    ax1.legend(loc="upper right")
    #ax1.set_title(params["ax1_title"])
    ax1.set_ylabel("error HRFE2 - HF")
    ax1.xaxis.set_ticks(p)
    ax1.set_xlim(params["xlimit"][0], params["xlimit"][1])
    ax1.yaxis.set_major_formatter(ticker.PercentFormatter(xmax=1, decimals=0))
    #ax4 = ax1.twinx()
    #ax4.yaxis.set_ticks([x for x in error_rom if x > params["ax2_ylimit"]])
    #ax4.yaxis.set_major_formatter(ticker.PercentFormatter(xmax=1, decimals=0))

    ax2.xaxis.set_ticks(p)
    ax2.set_xlim(params["xlimit"][0], params["xlimit"][1])
    ax2.set_ylim(0, params["ax2_ylimit"])
    ax2.yaxis.set_major_formatter(ticker.PercentFormatter(xmax=1, decimals=2))

    ax3.yaxis.set_ticks(error_rom)
    #ax3.set_ylabel("error ROM - HF")
    ax3.set_ylim(0, params["ax2_ylimit"])
    ax3.yaxis.set_major_formatter(ticker.PercentFormatter(xmax=1, decimals=2))


    plt.savefig(params["fig_fname"])
    plt.show()
