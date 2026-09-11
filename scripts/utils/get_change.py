"""dependency free script for getting percent change between round robin and similarity methods"""
import sys
def sim_change(x):
    res = ['0','0']
    for idx in range(2):
        num = (x[2]-x[idx])
        if num != 0:
            res[idx] = str((num/x[idx])*100)
    return ' '.join(res)

if __name__ == "__main__":
    pth = sys.argv[1]
    with open(pth, mode='r') as f:
        lines = f.readlines()
        sums =  [float(x.strip().split(' ')[1]) for x in lines[1:]]
        means =  [float(x.strip().split(' ')[2]) for x in lines[1:]]
        vars =  [float(x.strip().split(' ')[3]) for x in lines[1:]]
        print(sim_change(sums), sim_change(means), sim_change(vars))

