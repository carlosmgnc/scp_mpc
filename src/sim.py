import numpy as np
from sympy import *
from plot import Plots
from solve import solveProblem
from mpc import MPC

class simulation:
    def __init__(self, opt, solver):
        
        self.opt = opt
        self.solver = solver

        self.trajectory = np.zeros((14, self.opt.nk))
        self.trajectory[:, [0]] = self.solver.x[:, [0]]
        self.trajectory[1:4, [0]] = self.trajectory[1:4, [0]] + np.array([[0.5],[0],[0.5]])
        self.trajectory[4:7, [0]] = self.trajectory[4:7, [0]] + np.array([[0],[0],[0]])

        self.ol_trajectory = self.solver.x[:, [0]]

        self.u = self.solver.u
        self.sigma_opt = self.solver.sigma

        self.mpc = MPC(self.opt)
        self.u_mpc = np.zeros(3)

    # nonlinear dynamics derivative function for rk4
    def nonlin_func(self, t, X):
        x = X.flatten()

        alphak, betak, tk, tk1, k = self.opt.get_alphas(t)
        if k + 1 <= self.opt.nk-1:
            u = (alphak * self.u[:, k] + betak * self.u[:, k+1]).reshape((3, 1))
            self.u_mpc = (alphak * self.u1 + betak * self.u2)
            u1 = (alphak * self.u[:, k] + betak * self.u[:, k+1])
        else: 
            u = (self.u[:, k]).reshape((3, 1))
            self.u_mpc = self.u1
        u = u.flatten()


        def sigmaE_f(E_f):
            return self.sigma_opt * E_f

        return sigmaE_f(self.opt.E_f(x, u + self.u_mpc))
    
    # integrate nonlinear dynamics using rk4
    def integrate_full_trajectory(self):
        nsub = 30
        dt_sub = self.opt.dt / (nsub + 1)
        P_temp = self.trajectory[:, [0]]

        for i in range(0, self.opt.nk - 1):

            #self.mpc.canonicalize(self.trajectory[:, i], i)
            self.u1, self.u2 = self.mpc.solve_cvxpy(self.trajectory[:, i], i)
            print(i)

            print(self.u_mpc)
            for j in range(0, nsub + 1):
                sub_time = i * self.opt.dt + j * dt_sub
                P_temp = self.opt.rk41(self.nonlin_func, sub_time, P_temp, dt_sub)

            self.trajectory[:, [i + 1]] = P_temp

solver = solveProblem()
plotter = Plots()
solver.solve()

sim = simulation(solver.opt, solver)

sim.ol_trajectory = solver.x
sim.integrate_full_trajectory()

print("final position: " + str(sim.trajectory[1:4, -1]))
plotter.plot(sim.opt, solver, sim.trajectory)

