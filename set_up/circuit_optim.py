import qutip as qt 
import numpy as np
from scipy.optimize import minimize

sx     = qt.sigmax()
sy     = qt.sigmay()
sz     = qt.sigmaz()
si     = qt.identity(2)

class circuit: #generates circuit cells from http://cds.cern.ch/record/931003/files/0602174.pdf . n is the number of qubits
    def __init__(self, n, g, init_discrete = None, init_contin = None):
        # n: number of qubits
        # g: number of gadgets
        # du: derivatives
        # gamma: noise value
        # kraus: set of kraus operators 
        self.g = g
        self.n = n
        self.init_state = qt.Qobj(qt.basis(2**n,0), dims = [[2]*self.n, [1]*self.n ])

        #I think it might actually be worth setting gamma, du and kraus as global variables
        if init_discrete is None:
            self.sq  =  [np.random.choice(2) for i in range(g*n*(n+1))]
            self.tq  =  [np.random.choice(2) for i in range(g*n*(n-1))]
           
        else:
            self.sq, self.tq = init_discrete[0], init_discrete[1]

        if init_contin is None:
            self.sqp =  [2*np.pi*np.random.random() for i in range(int(3*sum(self.sq)) )]
        else:
             self.sqp = init_contin
        #size of sq and sqp is 2*g*n*(n+1)
        #size of tq is g*n*(n-1)
        print(self.sq, self.tq)
        self.cnots = []

        count = 0 
        for i in range(g):
            for j in range(n):
                temp_c = qt.tensor([si]*n)
                for k in range(n):
                    if k != j:
                        #print(self.tq[count])
                        if self.tq[count] == 1:
                            #print('cnot',qt.cnot(self.n, j, k))
                            temp_c = qt.cnot(self.n, j, k)*temp_c
                            count += 1
                            self.cnots.append(temp_c)
                        else:
                            self.cnots.append(temp_c)
                            count += 1

    def optimise(self):
        x0 = [np.random.random() for i in range(len(self.sqp))]
        res = minimize(circuit.func, x0, args=(self))
        self.sqp = res.x
        self.state = qt.Qobj(np.abs(self.unitary()*self.init_state), dims = [[2]*self.n, [1]*self.n ])

    def func(x0, *args):
        self = args[0]
        self.sqp = x0
        print(self.sqp)
        #print(self.unitary().dims,self.init_state.dims)
        self.state = qt.Qobj(np.abs(self.unitary()*self.init_state), dims = [[2]*self.n, [1]*self.n ])
        fid = np.real(1-np.abs(self.state.dag()*target_state))
        print(fid)
        return fid

    def unitary(self):
        count = 0
        c_count = 0
        p_c = 0
        #print(self.sqp) #need to fix the problem of empty array showing up here
        phases = np.array_split(self.sqp, sum(self.sq)) #fix this if sum(self.sq) = 0
        #print(phases)
        #print(len(self.sq),sum(self.sq), self.sq)
        #print(phases)
        u = qt.tensor([si]*self.n)
        for i in range(self.g):
            for j in range(self.n+1):
                squt = []
                for k in range(self.n):
                    if self.sq[count] == 1:
                        #print(phases[p_c])
                        h = phases[p_c][0]*sx + phases[p_c][1]*sy +phases[p_c][2]*sz
                        squt.append(qt.Qobj.expm(-1j*h))
                        p_c += 1
                    else:
                        squt.append(si)
                        #count += 1
                    count += 1
                    
                if j == self.n:
                    #print('final')
                    u = qt.tensor(squt)*u
                else:
                    #print(c_count)
                    #print(self.cnots[c_count],qt.tensor(squt))
                    u = self.cnots[c_count]*qt.tensor(squt)*u
                    c_count += 1
        self.uni = u
        return u


target_state = np.load('out_state.npy')
target_state = qt.Qobj(target_state, dims = [[2]*8,[1]*8]).unit()

n = 8
g = 1
sq = [0 for _ in range(n*(n+1))]

sq[0:8] = [1 for i in range(8)]

circ = circuit(n, g)
circ.optimise()

np.save('circ_full_State', circ.state)
np.save('one_qubit', circ.sq)
np.save('two_qubit', circ.tq)
np.save('phases', circ.sqp)

print(circ.sq)
print(circ.tq)
import matplotlib.pyplot as plt
fig, ax = qt.hinton(circ.state.ptrace((0,1,2,3)))
plt.show()














