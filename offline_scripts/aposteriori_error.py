import sys
import numpy as np
import argparse


def compute_mean_relative_error(comp, ref):
    length = np.shape(ref)[0]
    length_p = np.shape(comp)[0]
    if length_p < length:
            length = length_p
    ref = ref[:length, :]
    comp = comp[:length, :]
    n = np.shape(comp)[0] * np.ones((1, np.shape(comp)[1]))
    ERN1 = np.linalg.norm(comp / ref - ref / ref, ord=1, axis=0)
    return ERN1 / n


def compute_mean_normalized_error(comp, ref):
    length = np.shape(ref)[0]
    length_p = np.shape(comp)[0]
    if length_p < length:
            length = length_p
    ref = ref[:length, :]
    comp = comp[:length, :]
    EN1 = np.linalg.norm(comp - ref, ord=1, axis=0)
    N1 = np.linalg.norm(ref, ord=1)
    return EN1 / N1


def compute_mean_absolute_error(comp, ref):
    length = np.shape(ref)[0]
    length_p = np.shape(comp)[0]
    if length_p < length:
            length = length_p
    ref = ref[:length, :]
    comp = comp[:length, :]
    N = np.shape(comp)[0] * np.ones((1, np.shape(comp)[1]))
    return np.linalg.norm(comp - ref, ord=1, axis=0) / N


def compute_mean_square_error(comp, ref):
    length = np.shape(ref)[0]
    length_p = np.shape(comp)[0]
    if length_p < length:
            length = length_p
    ref = ref[:length, :]
    comp = comp[:length, :]
    N = np.shape(comp)[0] * np.ones((1, np.shape(ref)[1]))
    return np.linalg.norm(comp - ref, axis=0) / N


if __name__ == "__main__":
    #parser = argparse.ArgumentParser(description="Computes several error measurements of homogenized tension.")
    #parser.add_argument('ref_ht', help="reference homogenized tension filename (numpy format)")
    #parser.add_argument('comp_ht', help="computed homogenized tension filename (numpy format)")
    #args = parser.parse_args()
    if len(sys.argv) < 3:
        print("Computes several error measurements of homogenized tension.")
        print("arguments: reference measurement1 [measurement2 [measuremente3 [...]]]")
        sys.exit()
    compute_error = compute_mean_normalized_error
    #compute_error = compute_mean_square_error
    #compute_error = compute_mean_relative_error
    #compute_error = compute_mean_absolute_error
    np.set_printoptions(linewidth=1200, precision=3)
    f = sys.argv[1]
    ref = np.loadtxt(f)
    for f in sys.argv[2:]:
        print(f)
        measurement = np.loadtxt(f)
        print("{:.3e}".format(np.max(compute_error(measurement, ref)[6:])))

