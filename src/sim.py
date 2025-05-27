import numpy as np
from sympy import *
from plot import Plots
from solve import solveProblem
from mpc import MPC
import matplotlib.pyplot as plt
import time

class simulation:
    def __init__(self, solver, mpc):
        
        self.solver = solver
        self.mpc = mpc

        self.trajectory = np.zeros((14, self.solver.opt.nk))
        self.trajectory[:, [0]] = self.solver.x[:, [0]]

        self.ol_trajectory = self.solver.x[:, [0]]

        self.u = self.solver.u
        self.sigma_opt = self.solver.sigma

        self.u_mpc = np.zeros(3)
        self.u_rcs = 0
        self.u_mpc_list = np.zeros((3, self.solver.opt.nk))
        self.u_rcs_list = np.zeros(self.solver.opt.nk)


    # nonlinear dynamics derivative function for rk4
    def nonlin_func(self, t, X):
        x = X.flatten()

        alphak, betak, tk, tk1, k = self.solver.opt.get_alphas(t)
        if k + 1 <= self.solver.opt.nk-1:
            u = (alphak * self.u[:, k] + betak * self.u[:, k+1]).reshape((3, 1))
        else: 
            u = (self.u[:, k]).reshape((3, 1))
        
        u = u.flatten()

        def dyn(E_f):
            return self.sigma_opt * (E_f + self.solver.opt.B_rcs_ct * self.u_rcs)

        return dyn(self.solver.opt.E_f(x, u + self.u_mpc))
    
    # integrate nonlinear dynamics using rk4
    def simulate(self):
        nsub = 30
        dt_sub = self.solver.opt.dt / (nsub + 1)
        P_temp = self.trajectory[:, [0]]

        for i in range(0, self.solver.opt.nk - 1):
            print("k: " + str(i))

            #custom solver
            # t0 = time.perf_counter()
            # self.mpc.canonicalize(self.trajectory[:, i], i)
            # t1 = time.perf_counter()
            # self.u_mpc = self.mpc.solve_custom(i)
            # self.u_mpc_list[:, i] = self.u_mpc
            # print("solver time: " + str((t1-t0)*1000))

            # cvxpy solver
            t0 = time.perf_counter()
            self.u_mpc, self.u_rcs = self.mpc.solve_cvxpy(self.trajectory[:, i], i)
            t1 = time.perf_counter()
            ms = (t1-t0)*1000
            print("solver time: " + str(f"{ms:.2f}") + " ms")

            self.u_mpc_list[:, i] = self.u_mpc
            self.u_rcs_list[i] = self.u_rcs
            

            for j in range(0, nsub + 1):
                sub_time = i * self.solver.opt.dt + j * dt_sub
                P_temp = self.solver.opt.rk41(self.nonlin_func, sub_time, P_temp, dt_sub)

            self.trajectory[:, [i + 1]] = P_temp

solver = solveProblem()
plotter = Plots()
solver.solve()

mpc = MPC(solver.opt)
mpc.def_cvxpy_problem()

np.random.seed(10)

def rand_quat_perturbation(max_angle_rad):
    axis = np.random.randn(3)
    axis /= np.linalg.norm(axis)
    angle = np.random.choice([-max_angle_rad, max_angle_rad])
    w = np.cos(angle/2)
    xyz = axis * np.sin(angle/2)
    return np.concatenate(([w], xyz))

def quat_multiply(q1, q2):
    w1, x1, y1, z1 = q1
    w2, x2, y2, z2 = q2
    return np.array([
        w1*w2 - x1*x2 - y1*y2 - z1*z2,
        w1*x2 + x1*w2 + y1*z2 - z1*y2,
        w1*y2 - x1*z2 + y1*w2 + z1*x2,
        w1*z2 + x1*y2 - y1*x2 + z1*w2])

num_tests = 50
rand_pos_perturbation = 0.25*np.random.uniform(-1, 1, size=(3,num_tests))
rand_vel_perturbation = 0.1*np.random.uniform(-1, 1, size=(3,num_tests))
max_angle = np.deg2rad(10)

for i in range(num_tests):
    sim = simulation(solver, mpc)
    sim.trajectory[1:4, 0] += rand_pos_perturbation[:, i]
    sim.trajectory[4:7, 0] += rand_vel_perturbation[:, i]

    q_nom = sim.trajectory[7:11,0]
    q_delta = rand_quat_perturbation(max_angle)
    q_init = quat_multiply(q_delta, q_nom)
    q_init /= np.linalg.norm(q_init)
    sim.trajectory[7:11,0] = q_init

    mpc.init_params()

    sim.simulate()
    plotter.plot(solver, sim)

    print("final position: " + str(sim.trajectory[1:4, -1]))

plt.show(block=False)
plt.pause(1)
input()
plt.close()