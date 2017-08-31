import math
import numpy

def bateman(times,C0,lmbda):
    n = len(C0)
    out_arr = []
    for t in times:
        out = 0
        for i in range(n):
            ltot = 1.
            for j in range(i,n-1): ltot *= lmbda[j]
            sm = 0
            for j in range(i,n):
                denom = 1
                for p in range(i,n): 
                    if not p == j:
                        denom *= (lmbda[p]-lmbda[j])
                sm += math.exp(-lmbda[j]*t)/denom
            out += C0[i]*ltot*sm
        out_arr.append(out)

    return out_arr

if __name__=="__main__":

    thalf = numpy.array([0.00075, 0.001043379, 1.])
    lmbda = numpy.log(2)/thalf
    lmbda[2] = 0

    C0 = numpy.array([3.48e-4, 1e-30, 1e-30])
    times = numpy.arange(0.1,50.1,0.1)/365.

    CI = bateman(times,[C0[0]],[lmbda[0]])
    CXe = bateman(times,C0[0:2],lmbda[0:2])
    CCs = bateman(times,C0,lmbda)

        
