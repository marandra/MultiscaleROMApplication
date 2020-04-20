"""
Functions for data processsing
"""
import json
from pathlib import Path
import numpy
import pandas
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import aposteriori_error
from offline_common import Common


def plot_01():
    """
    Plot 1:
        Errors RVE case (ax1)
        Zoom errors RVE < 2% (ax2),
        ROM const error (ax3 inside ax2)
    """
    subset = ERRORS[(ERRORS["case"] == 5) & (ERRORS["points"] != "ROM")]
    errors = subset.pivot(index="points", columns="modes", values="error")
    errors_rom = ERRORS[(ERRORS["case"] == 5) & (ERRORS["points"] == "ROM")]

    # plot configuration
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
    fig.set_size_inches(PLOT_PARAMS["fig_size"][0], PLOT_PARAMS["fig_size"][1])

    for col in errors.columns:
        ax1.plot(errors.index, errors[col], label="{} modes".format(col), marker="")
        ax2.plot(errors.index, errors[col], label="{} modes".format(col), marker="o")

    for color, row in enumerate(errors_rom.itertuples(index=False)):
        ax3.plot(
            errors.index,
            [row.error] * len(errors.index),
            label=row.modes,
            marker="",
            linestyle="dotted",
            color="C{}".format(color),
        )

    ax1.set_title("Error HRFE$^2$ - HF")
    ax1.legend(loc="upper right")
    ax1.set_xticklabels("")
    ax1.xaxis.set_ticks(errors.index)
    ax1.set_xlim(PLOT_PARAMS["xlimit"][0], PLOT_PARAMS["xlimit"][1])
    ax1.yaxis.set_major_formatter(ticker.PercentFormatter(xmax=1, decimals=0))

    ax2.xaxis.set_ticks(errors.index)
    ax2.set_xlim(PLOT_PARAMS["xlimit"][0], PLOT_PARAMS["xlimit"][1])
    ax2.set_ylim(0, PLOT_PARAMS["ax2_ylimit"])
    ax2.yaxis.set_major_formatter(ticker.PercentFormatter(xmax=1, decimals=2))
    ax2.set_xlabel("Number of integration points")

    ax3.yaxis.set_ticks(errors_rom["error"])
    ax3.set_ylim(0, PLOT_PARAMS["ax2_ylimit"])
    ax3.yaxis.set_major_formatter(ticker.PercentFormatter(xmax=1, decimals=2))

    plt.savefig(
        PLOT_PARAMS["fig_fname"] + ".pdf", dpi=1000, bbox_inches="tight",
    )
    plt.show()


class Analysis(Common):
    """
    Functions for data processsing (mainly errors and times)
    """

    def __init__(self, **kargv):
        super(Analysis, self).__init__(**kargv)

    def compute_errors(self):
        """
        Returns DataFrame with (case, mode, points, error)
        """
        c_t = []
        c_m = []
        c_p = []
        c_e = []
        for t_id in self.context["cases_test_dataset"]:
            stress_ref = numpy.loadtxt(
                self.training_path / self.case_name(t_id) / "homogenized_stress.dat"
            )
            for m_id in self.context["rve_data_modes"]:
                for p_id in self.ip_subsets:
                    path = (
                        self.multiscale_path
                        / self.case_name(t_id)
                        / "_{}m_{}ip".format(m_id, p_id)
                    ).resolve()
                    stress = numpy.loadtxt(path / "homogenized_stress.dat")
                    err = max(
                        aposteriori_error.compute_mean_normalized_error(
                            stress, stress_ref
                        )[6:]
                    )
                    c_t.append(t_id)
                    c_m.append(m_id)
                    c_p.append(p_id)
                    c_e.append(err)
        return pandas.DataFrame(
            {"case": c_t, "modes": c_m, "points": c_p, "error": c_e}
        )

    # def read_time(t, ip, modes):
    #     with open("../training/validation/time_trajectory_{}/time.dat".format(t)) as fi:
    #         for line in fi.readlines():
    #             if "User time" in line:
    #                 utime = float(line.split(":")[-1])
    #                 continue
    #             if "System time" in line:
    #                 stime = float(line.split(":")[-1])
    #                 continue
    #         time_hf = utime + stime

    #     times = {}
    #     for m in modes:
    #         time_l = []
    #         for p in ip:
    #             path = "../multiscale_1ip_speedup/case_{}t_{}m_{}ip/time.dat".format(
    #                 t, m, p
    #             )
    #             with open(path) as fi:
    #                 for line in fi.readlines():
    #                     if "User time" in line:
    #                         utime = float(line.split(":")[-1])
    #                         continue
    #                     if "System time" in line:
    #                         stime = float(line.split(":")[-1])
    #                         continue
    #                 time_l.append(utime + stime)
    #         times[m] = time_l
    #     t_df = pandas.DataFrame(times, index=ip)
    #     return t_df, time_hf


##################################
# main
##################################
if __name__ == "__main__":

    import sys

    if len(sys.argv) > 1:
        DATA = Analysis(config_fname=sys.argv[1])
    else:
        DATA = Analysis()
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

    # plot configuration
    PLOT_PARAMS = json.loads(Path("./plot_params.json").read_text())
    # data preparation
    ERRORS = DATA.compute_errors()

    # Plot 1:
    #     Errors RVE case (ax1)
    #     Zoom errors RVE < 2% (ax2),
    #     ROM const error (ax3 inside ax2)
    plot_01()
