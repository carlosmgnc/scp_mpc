import numpy as np
import qp_module

class MPC:
    def __init__(self, opt):
        self.opt = opt
        self.N = 10

        self.n = self.N*(14 + 3)
        self.m = 14*(self.N)
        self.p = 7*(self.N - 1)

        self.Q = np.zeros((self.n, self.n))
        self.q = np.zeros(self.n)
        self.G = np.zeros((self.p, self.n))
        self.h = np.zeros(self.p)
        self.A = np.zeros((self.m, self.n))
        self.b = np.zeros(self.m)

        self.H = np.diag([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])
        self.F = np.diag([1, 1, 1])

    def canonicalize(self, x0, k):
        # form A matrix
        self.A[0:14, 0:(3+14)] = np.hstack([-self.opt.B.value[:, k*3: k*3 + 3], np.eye(14)])
        for i in range(0, self.N-1):
            indx_row = 14 + (i*14)
            indx_col = 3 + i*(14 + 3)

            indx_A = (k+i)*14
            indx_B = (k+i)*3

            A = self.opt.A.value[:, indx_A : indx_A + 14]
            B = self.opt.B.value[:, indx_B : indx_B + 3] + self.opt.C.value[:, indx_B : indx_B + 3]
            self.A[indx_row : indx_row + 14,  indx_col: indx_col + (14 + 3 + 14)] = np.hstack([-A, -B, -np.eye(14)])

        # form b vector
        self.b[0:14] = self.A[0:14, 0:14] @ x0
        
        #form G matrix
        tmax = self.opt.Tmax
        tmin = self.opt.Tmin
        
        for i in range(0, self.N-1):
            u = self.opt.u.value[:, [k + i]]
            fu_lower = np.linalg.norm(u) * u
            fu_lower = fu_lower.flatten()
        
            Fu = np.array([[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1], fu_lower])

            indx_row = i*Fu.shape[0]
            indx_col = i*(14 + 3)
            self.G[indx_row :indx_row + Fu.shape[0], indx_col : indx_col + Fu.shape[1]] = Fu

        # form h vector
        for i in range(0, self.N-1):
            indx_h = i*Fu.shape[0]
            h_i = np.array([tmax, tmax, tmax, tmax, tmax, tmax, tmin * (u.T@u)[0,0]])
            self.h[indx_h : indx_h + Fu.shape[0]] = h_i
        
        # form Q matrix
        for i in range(0, self.N):
            indx_Q = i*(14 + 3)
            self.Q[indx_Q: indx_Q + 14, indx_Q : indx_Q + 14] = self.H
            self.Q[indx_Q + 14 : indx_Q + 17, indx_Q + 14 : indx_Q + 17] = self.F
    
    def solve(self):
        print(self.Q.shape)
        print(self.q.shape)
        print(self.G.shape)
        print(self.h.shape)
        print(self.A.shape)
        print(self.b.shape)
        qp = qp_module.QP(self.Q, self.q, self.G, self.h, self.A, self.b)
        qp.solve()
        return qp.x

        