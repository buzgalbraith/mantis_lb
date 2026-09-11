"""dependincy free python code to get mean and variance"""
import sys
def mean(x):
    return sum(x)/len(x)
def var(x):
    m = mean(x)
    return sum((y - m)**2 for y in x) / len(x)
if __name__ == "__main__":
    pth = sys.argv[1]
    with open(pth, mode='r') as f:
        lines = f.readlines()
        lines = [int(x.strip()) for x in lines]
        print(mean(lines), var(lines))

