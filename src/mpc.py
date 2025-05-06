import numpy as np
import qp_module
import cvxpy as cvx

class MPC:
    def __init__(self, opt):
        self.opt = opt
        self.N = 35

        self.n = self.N*(14 + 3)
        self.m = 14*(self.N)
        self.p = 7*(self.N - 1)

        self.Q = np.zeros((self.n, self.n))
        self.q = np.zeros(self.n)
        self.G = np.zeros((self.p, self.n))
        self.h = np.zeros(self.p)
        self.A = np.zeros((self.m, self.n))
        self.b = np.zeros(self.m)

        self.H = np.diag([1, 100, 100, 100, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])
        self.H2 = np.diag([1, 10, 10, 10, 10, 10, 10, 1000, 1000, 1000, 1000, 1000, 1000, 1000])
        self.F = np.diag([1, 1, 1])

    def canonicalize(self, x0, k):

        self.A.fill(0)
        self.b.fill(0)
        # form A matrix
        B = self.opt.B.value[:, k*3: k*3 + 3] + self.opt.C.value[:, k*3: k*3 + 3]
        self.A[0:14, 0:(3+14)] = np.hstack([-B, np.eye(14)])
        
        for i in range(0, self.N-2):
            indx_row = 14 + (i*14)
            indx_col = 3 + i*(14 + 3)

            if k + i < self.N-2:
                indx_A = (k+i + 1)*14
                indx_B = (k+i + 1)*3
                A = self.opt.A.value[:, indx_A : indx_A + 14]
                B = self.opt.B.value[:, indx_B : indx_B + 3] + self.opt.C.value[:, indx_B : indx_B + 3]
                
                self.A[indx_row : indx_row + 14,  indx_col: indx_col + (14 + 3 + 14)] = np.hstack([-A, -B, np.eye(14)])
                
                x_refi = self.opt.x.value[:, k + i]
                urefi = self.opt.u.value[:, k + i]
                xref_i_1 =  self.opt.x.value[:, k + i + 1]
                self.b[i*14 + 14: i*14 + 14 + 14] = A@x_refi + B@urefi - xref_i_1
            
            else:
                self.A[indx_row : indx_row + 14,  indx_col: indx_col + (14 + 3 + 14)] = np.hstack([-np.eye(14), -np.zeros((14, 3)), np.eye(14)])
                self.b[indx_row : indx_row + 14] = np.zeros(14)


        # # form b vector
        indx_A = (k)*14
        indx_B = (k)*3
        A = self.opt.A.value[:, indx_A:indx_A + 14]
        B = self.opt.B.value[:, indx_B : indx_B + 3] + self.opt.C.value[:, indx_B : indx_B + 3]
        self.b[0:14] =  A @ x0 + B @ self.opt.u.value[:, k] - self.opt.x.value[:, k+1]

        
        # #form G matrix
        # tmax = self.opt.Tmax
        # tmin = self.opt.Tmin
        
        # for i in range(0, self.N-1):
        #     if k + i < self.opt.nk - 1:
        #         u = self.opt.u.value[:, [k + i]]
        #         fu_lower = np.linalg.norm(u) * u
        #         fu_lower = fu_lower.flatten()
            
        #         Fu = np.array([[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1], fu_lower])

        #         indx_row = i*Fu.shape[0]
        #         indx_col = i*(14 + 3)
        #         self.G[indx_row :indx_row + Fu.shape[0], indx_col : indx_col + Fu.shape[1]] = Fu

        # # form h vector
        # for i in range(0, self.N-1):
        #     indx_h = i*Fu.shape[0]
        #     h_i = np.array([tmax, tmax, tmax, tmax, tmax, tmax, tmin * (u.T@u)[0,0]])
        #     self.h[indx_h : indx_h + Fu.shape[0]] = h_i
        
        # form Q matrix
        for i in range(0, self.N):
            indx_Q = i*(14 + 3)
            self.Q[indx_Q: indx_Q + 3, indx_Q : indx_Q + 3] = self.F
            self.Q[indx_Q + 3 : indx_Q + 17, indx_Q + 3 : indx_Q + 17] = self.H
    
    def solve(self, k):
        qp = qp_module.QP(self.Q, self.q, self.G, self.h, self.A, self.b)
        qp.solve()
        print(qp.x[0:3])
        return  qp.x[0:3] + self.opt.u.value[:, k]
    
    def solve_cvxpy(self, x0, k):

        x = cvx.Variable((14, self.N))
        u = cvx.Variable((3, self.N))

        constraints = []
        cost = 0
        constraints += [x[:, 0] == x0[:] - self.opt.x.value[:, k]]

        for i in range(0, self.N-1):
            
            if k + i < self.N-1:
                cost += cvx.quad_form(x[:, i], self.H) + cvx.quad_form(u[:, i], self.F)
                indx_A = (k+i)*14
                indx_B = (k+i)*3
                A = self.opt.A.value[:, indx_A : indx_A + 14]
                B = self.opt.B.value[:, indx_B : indx_B + 3]
                C = self.opt.C.value[:, indx_B : indx_B + 3]

                constraints += [x[:, i+1] == A @(x[:, i]) + B @(u[:, i]) + C @(u[:, i+1])]

                uref_i = self.opt.u.value[:, k + i]
                constraints += [u[:, i] + uref_i <= self.opt.Tmax]
                constraints += [u[:, i] + uref_i >= -self.opt.Tmax]
                #constraints += [cvx.norm(u[:, i] + uref_i) <= self.opt.Tmax]
            else:
                cost += cvx.quad_form(x[:, i], self.H2) + cvx.quad_form(u[:, i], self.F)

        objective = cvx.Minimize(cost)
        prob = cvx.Problem(objective, constraints)
        prob.solve(solver="OSQP")

        return u.value[:, 0], u.value[:, 1]