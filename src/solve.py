import numpy as np
from scp import optProblem
import time

class solveProblem():
    def __init__(self):

        self.opt = optProblem()
        self.opt.def_cvx_problem()

        self.max_iter = 15
        self.converged_iter = self.max_iter
        self.sigma_list = np.empty((self.max_iter, 1))
        self.trajectory_list = np.zeros((self.max_iter, 14, self.opt.nk))

        self.x = np.zeros((14, self.opt.nk))
        self.u = np.zeros((3, self.opt.nk))
        self.sigma = 0
        self.cost = 0
        self.nu = np.zeros((14, self.opt.nk-1))

    # sequentially solve convex sub-problems
    def solve(self):
        for i in range(0, self.max_iter):

            t0 = time.perf_counter()
            self.opt.discretize()
            t1 = time.perf_counter()
            ms = (t1-t0)*1000
            print("discretize time: " + str(f"{ms:.2f}") + " ms")

            self.x, self.u, self.sigma, self.cost, self.nu = self.opt.solve_cvx_problem()

            self.sigma_list[i] = self.sigma
            self.trajectory_list[i, :, :] = self.x

            if np.abs(self.sigma - self.cost) <= 0.00025:
                self.converged_iter = i
                break
        return 