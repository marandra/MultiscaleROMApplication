import os
import time
import main_kratos_hprom as kratos

TRAJ = []
TRAJ.append([0.001, 0.0,   0.0, ])  # 0
TRAJ.append([0.0,   0.001, 0.0, ])  # 1
TRAJ.append([0.0,   0.0,   0.001])  # 2
TRAJ.append([0.001, 0.001, 0.0, ])  # 3
TRAJ.append([0.001, 0.0,   0.001])  # 4
TRAJ.append([0.0,   0.001, 0.001])  # 5
TRAJ.append([0.001, 0.001, 0.001])  # 6


if __name__ == "__main__":

    if False:
        # HPROM
        for m in [10]:
            for p in [100]:
                filename = "bases/rve_{}_{}.json".format(m, p)
                k = kratos.Kratos(filename=filename)
                #for t in range(63):
                for t in  [1]:
                    path_out ="traj_{}/m{}p{}".format(t, m, p)
                    if not os.path.exists(path_out):
                        os.makedirs(path_out)
                    t1 = time.time()
                    k.init_case(strain=TRAJ[t], path_out=path_out)
                    k.run()                                     
                    t2 = time.time()
                    with open("{}/time.dat".format(path_out), 'w') as fo:
                        fo.write("{:.2f}\n".format(t2 - t1))

    if True:
        # ROM
        for m in [10]:
            for p in ['ROM']:
                filename = "bases/rve_{}_{}.json".format(m, p)
                k = kratos.Kratos(filename=filename)
                for t in [1]:
                    path_out ="traj_{}/m{}p{}".format(t, m, p)
                    if not os.path.exists(path_out):
                        os.makedirs(path_out)
                    t1 = time.time()
                    k.init_case(strain=TRAJ[t], path_out=path_out)
                    k.run()                                     
                    t2 = time.time()
                    with open("{}/time.dat".format(path_out), 'w') as fo:
                        fo.write("{:.2f}\n".format(t2 - t1))
