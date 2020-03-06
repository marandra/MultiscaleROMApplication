import aposteriori_error
import numpy
from pathlib import Path, PurePath
import glob
import pprint
import matplotlib.pyplot as plt


def read_data(stress_ref):
    # base_path = Path("../multiscale_1ip")
    case_paths = glob.glob("../multiscale_1ip/case_*/homogenized_stress.dat")
    errors = {}
    for case_path in case_paths:
        stress = numpy.loadtxt(case_path)
        case = case_path.split("/")[2]
        errors[case] = max(
            aposteriori_error.compute_mean_normalized_error(stress, stress_ref)[6:]
        )
    for c in sorted(errors, key=errors.get):
        print(c, errors[c])
        # pprint.pprint(kc, errors[kc])
    return errors


def make_matrix(errors):
    m = []
    p = []
    for k, v in errors.items():
        mi = int(k.split("_")[2].split("m")[0])
        m.append(mi)
        pi = int(k.split("_")[3].split("ip")[0])
        p.append(pi)
    print(m)
    print(p)
    ms = list(sorted(set(m)))
    ps = list(sorted(set(p)))
    a = {}
    print(ms)
    print(ps)
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
    stress_ref = numpy.loadtxt(
        "../training/validation/trajectory_39/homogenized_stress.dat"
    )
    errors = read_data(stress_ref)
    # errors_series = make_matrix(errors)
    a, p = make_matrix(errors)
    print(a)
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=False)
    for k, v in a.items():
        ax1.plot(p, v, label=k, marker="")
        ax2.plot(p, v, label=k, marker="o")
    # configure plot
    ax2.legend()
    ax2.set_ylim(0, 0.0001)
    # ax1.spines['bottom'].set_visible(False)
    # ax2.spines['top'].set_visible(False)
    # ax1.xaxis.tick_top()
    # ax1.tick_params(labeltop=False)  # don't put tick labels at the top
    # ax2.xaxis.tick_bottom()

    plt.show()
