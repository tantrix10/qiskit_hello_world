import numpy as np 
from scipy.optimize import minimize as m 
import matplotlib.pyplot as plt 
import qutip as qt
np.set_printoptions(linewidth= np.inf)

goal = np.load('hello_world.npy')
goal_tri = np.triu(goal,1)

#import the goal state and just make sure it looks like it should!
# print(goal)

# plt.imshow(goal_tri)
# plt.show()

def objective(x):
	state = np.triu(state_gen(x).ptrace((0,1,2,3)).full(),1)
	return m(norm_norm, np.random.random(), args = (state)).fun
	#return np.linalg.norm(state- goal)


def norm_norm(x,*args):
	return np.linalg.norm(args[0]-x*goal_tri)

def state_gen(x):
	#state = np.random.random(16)
	state = qt.Qobj(x, dims = [[2]*n, [1]*n])
	#print(state,np.dot(state,state))
	return state


# print(objective(np.random.random(16)))

n = 8

res = m(objective, np.random.random(2**n))

print(res.fun)
print(np.linalg.norm(goal_tri))

out_state = state_gen(res.x)


np.save('out_state', np.array(res.x))

# np.save('out_state', out_state)
# outer = np.outer(out_state,out_state)

# print(out_state)
# print(outer)
# print(np.triu(outer,1))

plt.imshow(np.real(np.triu(state_gen(out_state).ptrace((0,1,2,3)).full(),1)))
plt.show()