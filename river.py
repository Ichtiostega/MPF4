import pickle as pckl
import numpy as np 
from copy import deepcopy
from matplotlib import pyplot as plt 

class River:
    def gen_mat(self):
        self.AA = np.zeros((self.size,self.size))
        self.BB = np.zeros((self.size,self.size))

        for i in range(self.size-1):
            self.AA[i][i] = 1 + self.C_d
            self.BB[i][i] = 1 - self.C_d
            self.AA[i][i+1] = self.C_a/4 - self.C_d/2 
            self.BB[i][i+1] = self.C_d/2 - self.C_a/4
            self.AA[i+1][i] = - self.C_a/4 - self.C_d/2 
            self.BB[i+1][i] = self.C_d/2 + self.C_a/4

        self.AA_inv = np.linalg.inv(self.AA)
        self.AA_inv_BB = np.matmul(self.AA_inv, self.BB)


    def __init__(self, dt=0.1, dx=1, n=10, size=1000, U=0.1, D=1):
        self.size = int(size * dx)
        self.body = np.zeros(self.size, dtype=float)
        self.l_in = int(10 * dx)
        self.l_out = int(250 * dx)
        self.A = 1.5
        self.U = U
        self.D = D
        self.m = 10000

        self.dt = dt
        self.dx = dx
        self.n = n

        self.C_in = self.m / (self.A * self.U * self.dt * self.n)
        self.C_a = self.U * self.dt / self.dx
        self.C_d = self.D * self.dt / self.dx**2
        print(self.C_a, self.C_d)
        self.gen_mat()

        self.iter = 0

        self.q_1 = ( self.C_d * (1 - self.C_a) - (self.C_a/6) * ( self.C_a**2 - 3*self.C_a + 2 ) )
        self.q_2 = ( self.C_d * (2 - 3*self.C_a) - (self.C_a/2) * ( self.C_a**2 - 2*self.C_a - 1 ) )
        self.q_3 = ( self.C_d * (1 - 3*self.C_a) - (self.C_a/2) * ( self.C_a**2 - self.C_a - 2 ) )
        self.q_4 = ( self.C_d * self.C_a + (self.C_a/6) * ( self.C_a**2 - 1 ) )
        print(self.q_1, self.q_2, self.q_3, self.q_4)

    def __str__(self):
        return str(self.body)

    def plot(self):
        plt.plot(list(range(len(self.body))), self.body)
        plt.show()

    def iterate_quickest(self):
        def _it():
            tmp = deepcopy(self.body)
            for i in range(2,len(self.body)-1):
                tmp[i] =    self.body[i] + self.body[i+1] * self.q_1 - self.body[i] * self.q_2 + \
                            self.body[i-1] * self.q_3 + self.body[i-2] * self.q_4
            return tmp

        if self.iter < self.n:
            print(self.iter, self.n, self.C_in)
            self.body[self.l_in] += self.C_in

        tmp = _it()

        while tmp[-1] > 0.00001:
            self.size += 100
            self.body = np.append(self.body, np.zeros(100, dtype=float))
            tmp = _it()
        
        self.iter += 1
        self.body = tmp

    def iterate_CN(self):
        def _it():
            tmp = deepcopy(self.body)
            tmp = np.matmul(self.AA_inv_BB, self.body)
            return tmp

        self.body[self.l_in] = self.body[self.l_in] + self.C_in if self.iter < self.n else self.body[self.l_in]

        tmp = _it()

        print(tmp[-1], tmp[-2], tmp[-3])
        while tmp[-3] > 0.00001:
            self.size += 100
            self.gen_mat()
            self.body = np.append(self.body, np.zeros(100, dtype=float))
            tmp = _it()
        
        self.iter += 1
        self.body = tmp


# r = River(dt=0.1)
# print(r.C_a, r.C_d)
# for i in range(100000):
#     r.iterate_quickest()
#     if i%10000 == 0:
#         print(sum(r.body))
#         r.plot()

# print(r)
# print(sum(r.body))
# r.plot()

# Zad 2

ns = [1, 100000]
out = [[n,{"Q": [], "CN": []}] for n in ns]
print(out)
for it in out:
    print(it)
    Q = River(dt=0.1, n=it[0])
    CN = River(dt=0.1, n=it[0])
    print(f'N: {it[0]}')
    i = 0
    fin = False
    reached_lout = False
    while not fin:
        Q.iterate_quickest()
        CN.iterate_CN()
        it[1]["Q"].append(Q.body[Q.l_out])
        it[1]["CN"].append(CN.body[CN.l_out])
        if not reached_lout:
            if it[1]["Q"][-1] > 0.0001 and it[1]["CN"][-1] > 0.0001:
                print("LOUT")
                reached_lout = True
        else:
            if it[1]["Q"][-1] < 0.0001 and it[1]["CN"][-1] < 0.0001:
                fin = True
        i += 1
        if i%1000 == 0:
            print(it[1]["Q"][-1], it[1]["CN"][-1])

with open("data.pckl", "wb") as f:
    pckl.dump(out, f)

# out = []
# with open("data.pckl", "rb") as f:
#     out = pckl.load(f)

# for it in out:
#     Qout = it[1]["Q"]
#     CNout = it[1]["CN"]
#     n = it[0]
#     print(n)
#     plt.plot(range(len(Qout)), Qout)
#     plt.show()
#     plt.plot(range(len(CNout)), CNout)
#     plt.show()

# zad 3

Us = [3, 2] # 2 - macierz jednostkowa
U_out = [[U, []] for U in Us]
for U in U_out:
    try:
        r = River(U=U[0])
        r.iterate_CN()
        U[1].append(deepcopy(r.body))
        r.plot()
        for i in range(5000):
            r.iterate_CN()
        U[1].append(deepcopy(r.body))
        r.plot()
    except Exception as e:
        print(e)

with open("data_U.pckl", "wb") as f:
    pckl.dump(U_out, f)

# U_out = []
# with open("data_U.pckl", "rb") as f:
#     U_out = pckl.load(f)

# for U in U_out:
#     print(U[0])
#     for data in U[1]:
#         plt.plot(range(len(data)), data)
#         plt.show()

Ds = [16, 0.05] # 16

D_out = [[D, []] for D in Ds]
for D in D_out:
    try:
        r = River(D=D[0])
        r.iterate_CN()
        D[1].append(deepcopy(r.body))
        r.plot()
        for i in range(50000):
            r.iterate_CN()
        D[1].append(deepcopy(r.body))
        r.plot()
    except Exception as e:
        print(e)

with open("data_D.pckl", "wb") as f:
    pckl.dump(D_out, f)

# D_out = []
# with open("data_D.pckl", "rb") as f:
#     D_out = pckl.load(f)

# for D in D_out:
#     print(D[0])
#     for data in D[1]:
#         plt.plot(range(len(data)), data)
#         plt.show()