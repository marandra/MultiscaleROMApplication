import sys

with open(sys.argv[1], "r") as fi:
    i = 0
    for line in fi.readlines():
        if line[0] == "#":
            continue
        rstrain = line.split()[:6]
        l = ""
        for e, s in enumerate(rstrain):
            l += "rstrain({:>3},{})={}; ".format(i, e, s)
        print(l)
        i += 1
    l += "\n"

with open(sys.argv[1], "r") as fi:
    i = 0
    for line in fi.readlines():
        if line[0] == "#":
            continue
        rstress = line.split()[6:]
        l = ""
        for e, s in enumerate(rstress):
            l += "rstress({:>3},{})={}; ".format(i, e, s)
        print(l)
        i += 1
    l += "\n"

with open(sys.argv[2], "r") as fi:
    i = 0
    for line in fi.readlines():
        if line[0] == "#":
            continue
        rcm = line.split()
        l = ""
        for e, s in enumerate(rcm):
            l += "rcm({:>3},{:>2})={}; ".format(i, e, s)
            if e % 6 == 5:
                l += "\n"
        print(l)
        i += 1
