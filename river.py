import numpy as np 
from copy import deepcopy
from matplotlib import pyplot as plt 

class River:
    def __init__(self, dt=0.01, dx=1, n=10):
        self.body = np.zeros(1000, dtype=float)
        self.l_in = 10
        self.l_out = 250
        self.A = 1.5
        self.U = 0.1
        self.D = 1
        self.m = 10

        self.dt = dt
        self.dx = dx
        self.n = n

        self.C_in = self.m / (self.A * self.U * self.dx * self.n)
        self.C_a = self.U * self.dt / self.dx
        self.C_d = self.D * self.dt / self.dx**2

        self.iter = 0

    def __str__(self):
        return str(self.body)

    def iterate(self):
        self.body[self.l_in] = self.C_in if self.iter < self.n else 0
        tmp = deepcopy(self.body)
        def _it(tmp):
            for i in range(2,len(self.body)-1):
                tmp[i] =    self.body[i] + \
                            self.body[i+1] * ( self.C_d * (1 - self.C_a) - (self.C_a/6) * ( self.C_a**2 - 3*self.C_a + 2 ) ) + \
                            self.body[i] * ( self.C_d * (2 - 3*self.C_a) - (self.C_a/2) * ( self.C_a**2 - 2*self.C_a - 1 ) ) + \
                            self.body[i-1] * ( self.C_d * (1 - 3*self.C_a) - (self.C_a/2) * ( self.C_a**2 - self.C_a - 2 ) ) + \
                            self.body[i-2] * ( self.C_d * self.C_a + (self.C_a/6) * ( self.C_a**2 - 1 ) )
            return tmp

        tmp = _it(tmp)

        while tmp[-1] != 0:
            np.append(self.body, np.zeros(100, dtype=float))
            tmp = deepcopy(self.body)
            tmp[self.l_in] = self.C_in if self.iter <= self.n else 0
            tmp = _it(tmp)
        
        self.iter += 1
        self.body = tmp


r = River(dt=0.1)
for i in range(100000):
    r.iterate()
    if i%1000 == 0:
        print(r)
        plt.plot(list(range(len(r.body))), r.body)
        plt.show()

print(r)

plt.plot(list(range(len(r.body))), r.body)
plt.show()