import numpy as np

class MPC:
    def __init__(self, opt):
        self.opt = opt
        self.N = 10

        self.n = self.N*(14 + 3) - 1
        self.m = 14*self.N
        self.p = 1

        self.Q = np.zeros((self.n, self.n))
        self.q = np.zeros(self.n)
        self.G = np.zeros((self.p, self.n))
        self.h = np.zeros(self.p)
        self.A = np.zeros((self.m, self.n))
        self.b = np.zeros(self.n)
        # self.x0 = np.zeros(14)

    def canonicalize(self, k):
        # form A matrix
        self.A[0:14, 0:(3+14)] = np.hstack([-self.opt.B[k], np.eye(14)])
        for i in range(0, self.n):
            indx_row = 14 + (i*14)
            indx_col = 3 + i*(14 + 3)
            self.A[indx_row : indx_row + 14,  indx_col: indx_col + (14 + 3) ] = np.hstack[[-self.opt.A[k + i], -self.opt.B[k + i], -np.eye(14)]]
        
        #form G matrix
        tmax = self.opt.Tmax
        fu_lower = np.linalg.norm(self.opt.u[:, [k]]) * self.opt.u[:, [k]].T
        Fu = np.array([[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1], fu_lower])
        for i in range(0, self.n):
            indx_row = (i*Fu.shape[0])
            indx_col = i*(14 + 3)
            self.G[indx_row :indx_row + Fu.shape[0], indx_col : indx_col + Fu.shape[1]] = Fu
            
                

        
