import numpy as np
import time
from datetime import timedelta
import matplotlib.pyplot as plt
from variables import *

def make_random_herm(size):
    t1=np.random.rand(size,size)
    t2=1j*np.random.rand(size,size)
    test=np.zeros((size,size),dtype=complex)
    test=np.asmatrix(t1+t2)+np.asmatrix(t1+t2).conj().T
    return test

sizes=[]
diagtimes=[]
for i in np.linspace(10000,10000,1):
    t=make_random_herm(size=int(i))
    start_time=time.monotonic()
    np.linalg.eigh(t)
    end_time=time.monotonic()
    run_time=timedelta(end_time-start_time)
    print(f"{i/100}"+r'% done')
    diagtimes.append(run_time.seconds)
    sizes.append(i)

print(sizes,diagtimes)
# plt.scatter(sizes,diagtimes)
# plt.scatter(np.log(sizes),np.log(diagtimes))

# plt.show()
