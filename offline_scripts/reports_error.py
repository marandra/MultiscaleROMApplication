#from pathlib import Path, PurePath
import pprint
import json
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy
import pandas
import aposteriori_error


def read_time(t, ip, modes):
    with open("../training/validation/time_trajectory_{}/time.dat".format(t)) as fi:
        for line in fi.readlines():
            if "User time" in line:
                utime = float(line.split(":")[-1])
                continue
            if "System time" in line:
                stime = float(line.split(":")[-1])
                continue
        time_hf = utime + stime

    times = {}
    for m in modes:
        time_l = []
        for p in ip:
            path = "../multiscale_1ip_speedup/case_{}t_{}m_{}ip/time.dat".format(
                t, m, p
            )
            with open(path) as fi:
                for line in fi.readlines():
                    if "User time" in line:
                        utime = float(line.split(":")[-1])
                        continue
                    if "System time" in line:
                        stime = float(line.split(":")[-1])
                        continue
                time_l.append(utime + stime)
        times[m] = time_l
    t_df = pandas.DataFrame(times, index=ip)
    return t_df, time_hf


def read_data(t, ip, modes):
    stress_ref = numpy.loadtxt(
        "../training/skip_trajectory_{}/homogenized_stress.dat".format(t)
    )
    errors = {}
    for m in modes:
        em = []
        for p in ip:
            path = "../multiscale_1ip/case_{}t_{}m_{}ip/homogenized_stress.dat".format(
                t, m, p
            )
            stress = numpy.loadtxt(path)
            e = max(
                aposteriori_error.compute_mean_normalized_error(stress, stress_ref)[6:]
            )
            em.append(e)
        errors[m] = em
    e_df = pandas.DataFrame(errors, index=ip)
    return e_df


def read_rom(t, modes):
    stress_ref = numpy.loadtxt(
        "../training/skip_trajectory_{}/homogenized_stress.dat".format(t)
    )
    errors = []
    for m in modes:
        path = "../multiscale_1ip/ROM_cases/case_{}t_{}m_{}ip/homogenized_stress.dat".format(
            t, m, "ROM"
        )
        stress = numpy.loadtxt(path)
        e = max(aposteriori_error.compute_mean_normalized_error(stress, stress_ref)[6:])
        errors.append(e)
    e_s = pandas.Series(errors, index=modes)
    return e_s


##################################
# main
##################################
if __name__ == "__main__":

    # load parameters
    fp = open("config.json", "r")
    params = json.load(fp)
    pprint.pprint(params)
    fo = open(params["dat_fname"], "w")

    param_general = json.load(open("../configuration.json"))
    ip = param_general["config_data"]["rve_data_points"][:-1]  # skip ROM
    modes = param_general["config_data"]["rve_data_modes"]
    trajectory = param_general["config_data"]["skip_cases"][0]

    # compute errors HRFE2
    errors = read_data(trajectory, ip, modes)
    pprint.pprint(errors)
    print()

    # compute errors ROM
    errors_rom = read_rom(trajectory, modes)
    pprint.pprint(errors_rom)
    print()

    ## compute times
    # times, time_hf = read_time(trajectory, ip, modes)
    # pprint.pprint(times)
    # times.plot()
    # print()

    ## compute speedup
    # speedup = pandas.DataFrame(time_hf / times)
    # pprint.pprint(speedup)
    # speedup.plot()
    # print()

    # line = "Time HF: {}s".format(time_hf)
    # print(line)
    # fo.write(line + "\n")
    # first = True
    # for c in sorted(errors, key=errors.get):
    #    if first:
    #        er = errors[c]
    #        tr = times[c]
    #        first = False
    #        line = "{:<25} {}      {}      {}  {}  {}  {}".format(
    #            "case", "error", "time", "s/err", "s=tr/t", "err=e/er", "speedup"
    #        )
    #        print(line)
    #        fo.write(line + "\n")
    #        line = "----------------------------------------------------------------"
    #        print(line)
    #        fo.write(line + "\n")
    #    e = errors[c]
    #    t = times[c]
    #    s = tr / t
    #    err = e / er
    #    line = "{:<25} {:>1.3e} {:>6.2f}s {:>8.3f} {:>8.3f} {:>8.3f} {:>8.3f}".format(
    #        c, e, t, s / err, s, err, time_hf / t
    #    )
    #    print(line)
    #    fo.write(line + "\n")

    # configure plot
    # matplotlib.style.use("seaborn-colorblind")
    # matplotlib.style.use("seaborn")
    # matplotlib.style.use("bmh")
    plt.rcParams["text.latex.preamble"] = [r"\usepackage{lmodern}"]
    plt_params = {
        "text.usetex": True,
        "font.size": 10,
        "font.family": "serif",
    }
    plt.rcParams.update(plt_params)

    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=False)
    ax3 = ax2.twinx()
    fig.set_size_inches(params["fig_size"][0], params["fig_size"][1])

    for col in errors:
        ax1.plot(errors.index, errors[col], label="{} modes".format(col), marker="")
        ax2.plot(errors.index, errors[col], label="{} modes".format(col), marker="o")

    for color, e in enumerate(errors_rom):
        ax3.plot(
            errors.index,
            [e] * len(errors.index),
            label=e,
            marker="",
            linestyle="dotted",
            color="C{}".format(color),
        )

    ax1.set_title("Error HRFE$^2$ - HF")
    ax1.legend(loc="upper right")
    ax1.set_xticklabels("")
    ax1.xaxis.set_ticks(errors.index)
    ax1.set_xlim(params["xlimit"][0], params["xlimit"][1])
    ax1.yaxis.set_major_formatter(ticker.PercentFormatter(xmax=1, decimals=0))

    ax2.xaxis.set_ticks(errors.index)
    ax2.set_xlim(params["xlimit"][0], params["xlimit"][1])
    ax2.set_ylim(0, params["ax2_ylimit"])
    ax2.yaxis.set_major_formatter(ticker.PercentFormatter(xmax=1, decimals=2))
    ax2.set_xlabel("Number of integration points")

    ax3.yaxis.set_ticks(errors_rom.values)
    ax3.set_ylim(0, params["ax2_ylimit"])
    ax3.yaxis.set_major_formatter(ticker.PercentFormatter(xmax=1, decimals=2))

    plt.savefig(
        params["fig_fname"] + ".pdf", dpi=1000, bbox_inches="tight",
    )
    plt.show()
