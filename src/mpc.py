import numpy as np
import qp_module
import cvxpy as cvx
from cvxpygen import cpg
from mpc_socp_solver.cpg_solver import cpg_solve

class MPC:
    def __init__(self, opt):
        self.opt = opt
        self.N = 35

        self.n = (self.N - 1)*(14 + 3)
        self.m = 14*(self.N - 1)
        self.p = 6*(self.N - 1)

        self.H_vec = np.array([1, 1000, 1000, 1000, 1, 1, 1, 500, 500, 500, 500, 1, 1, 1])
        self.H2_vec = np.array([1, 10000, 10000, 10000, 1, 1, 1, 1000, 1000, 1000, 1000, 10, 10, 10])

        self.H = np.diag(self.H_vec)
        self.H2 = np.diag(self.H2_vec)
        self.F = np.diag([0, 0, 0])

        self.x = cvx.Variable((14, self.N))
        self.u = cvx.Variable((3, self.N-1))
        self.u_rcs = cvx.Variable(self.N-1)

        self.x0 = cvx.Parameter(14)

        self.x_ref = cvx.Parameter((14, self.N))
        self.u_ref = cvx.Parameter((3, self.N-1))

        self.A = cvx.Parameter((14, 14 * (self.N - 1)))
        self.B = cvx.Parameter((14, 3 * (self.N - 1)))
        self.B_rcs = cvx.Parameter((14, (self.N - 1)))
        self.H_sqrt = cvx.Parameter((14, self.N-1), nonneg=True)

        self.Q_custom = np.zeros((self.n, self.n))
        self.q_custom = np.zeros(self.n)
        self.G_custom = np.zeros((self.p, self.n))
        self.h_custom = np.zeros(self.p)
        self.A_custom = np.zeros((self.m, self.n))
        self.b_custom = np.zeros(self.m) 
    
    def def_cvxpy_problem(self):
        constraints = []
        cost = 0

        constraints += [self.x[:, 0] == self.x0 - self.x_ref[:, 0]]
        for i in range(0, self.N-1):
            
            indx_A = (i)*14
            indx_B = (i)*3
            A = self.A[:, indx_A : indx_A + 14]
            B = self.B[:, indx_B : indx_B + 3]

            constraints += [self.x[:, i+1] == A @ (self.x[:, i]) + B @ (self.u[:, i]) + self.B_rcs[:, i]*self.u_rcs[i]]
            constraints += [self.u_rcs[i] <= 0.5]
            constraints += [self.u_rcs[i] >= -0.5]

            constraints += [cvx.norm(self.u[:, i] + self.u_ref[:, i]) <= self.opt.Tmax]
            constraints += [self.opt.cos_delta_max.value * cvx.norm(self.u_ref[:, i] + self.u[:, i]) <= self.u_ref[0, i] + self.u[0, i]]

            cost += cvx.sum_squares(cvx.diag(self.H_sqrt[:, i]) @ self.x[:, i+1]) + cvx.quad_form(self.u[:, i], self.F)
        
        objective = cvx.Minimize(cost)
        self.prob = cvx.Problem(objective, constraints)

        #cpg.generate_code(self.prob, code_dir='mpc_socp_solver', solver='QOCO', prefix="mpc_socp_")
        self.prob.register_solve("mpc_socp", cpg_solve)


    def solve_cvxpy(self, x0, k):
        self.x0.value = x0

        #self.prob.solve(solver="QOCO")
        self.prob.solve(method="mpc_socp")
        print(self.prob.status)

        self.update_params()

        return self.u.value[:, 0], self.u_rcs.value[0]
    
    def init_params(self):
        self.x_ref.value = self.opt.x.value
        self.u_ref.value = self.opt.u.value[:, :-1]

        self.A.value = self.opt.A.value
        self.B_zoh = self.opt.B.value + self.opt.C.value
        self.B.value = self.B_zoh
        self.B_rcs.value = self.opt.B_rcs
        
        self.H_sqrt.value = np.tile(np.sqrt(self.H_vec).reshape(14, 1), (1, self.N-1))
        self.H_sqrt.value[:, -1] = np.sqrt(self.H2_vec)

    def update_params(self):
        new_x_ref = np.zeros((14, self.N))
        new_x_ref[:, :-1] = self.x_ref.value[:, 1:]
        new_x_ref[:, -1] = self.x_ref.value[:, -2]
        self.x_ref.value = new_x_ref
        
        new_u_ref = np.zeros((3, self.N-1))
        new_u_ref[:, :-1] = self.u_ref.value[:, 1:]
        new_u_ref[:, -1] = 0
        self.u_ref.value = new_u_ref
        
        new_A = np.zeros((14, 14*(self.N-1)))
        new_A[:, :-14] = self.A.value[:, 14:]
        new_A[:, -14:] = np.eye(14)
        self.A.value = new_A

        new_B = np.zeros((14, 3*(self.N-1)))
        new_B[:, :-3] = self.B.value[:, 3:]
        new_B[:, -3:] = 0
        self.B.value = new_B

        new_B_rcs = np.zeros((14, self.N-1))
        new_B_rcs[:, :-1] = self.B_rcs.value[:, 1:]
        new_B_rcs[:, -1] = 0
        self.B_rcs.value = new_B_rcs
        
        new_H = np.zeros((14, self.N-1))
        new_H[:, :-1] = self.H_sqrt.value[:, 1:]
        new_H[:, -1] = 0
        self.H_sqrt.value = new_H

    def canonicalize(self, x0, k):

        self.A.fill(0)
        self.b.fill(0)
        self.G.fill(0)
        self.h.fill(0)
        # form A matrix
        
        #initial condition dynamics constraint
        B = self.opt.B.value[:, k*3: k*3 + 3] + self.opt.C.value[:, k*3: k*3 + 3]
        self.A[0:14, 0:(3+14)] = np.hstack([-B, np.eye(14)])
        
        # remaining dynamics constraints
        for i in range(1, self.N-1):
            indx_row = (i*14)
            indx_col = 3 + (i-1)*(14 + 3)

            if k + i < self.N-1:
                indx_A = (k+i)*14
                indx_B = (k+i)*3
                A = self.opt.A.value[:, indx_A : indx_A + 14]
                B = self.opt.B.value[:, indx_B : indx_B + 3] + self.opt.C.value[:, indx_B : indx_B + 3]
                
                self.A[indx_row : indx_row + 14,  indx_col: indx_col + (14 + 3 + 14)] = np.hstack([-A, -B, np.eye(14)])
            
            else:
                self.A[indx_row : indx_row + 14,  indx_col: indx_col + (14 + 3 + 14)] = np.hstack([-np.eye(14), -np.zeros((14, 3)), np.eye(14)])

        # form b vector
        indx_A = (k)*14
        indx_B = (k)*3
        A = self.opt.A.value[:, indx_A:indx_A + 14]
        self.b[0:14] = A @ (x0[:] - self.opt.x.value[:, k]) 
        
        #form G matrix
        tmax = self.opt.Tmax
        tmin = self.opt.Tmin
        
        for i in range(0, self.N-1):
            if k + i < self.opt.nk - 1:
                u = self.opt.u.value[:, [k + i]]
                # fu_lower = np.linalg.norm(u) * u
                # fu_lower = fu_lower.flatten()
                fu_lower = np.array([0, 0, 0])
            
                Fu = np.array([[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]])

                indx_row = i*Fu.shape[0]
                indx_col = i*(14 + 3)
                self.G[indx_row :indx_row + Fu.shape[0], indx_col : indx_col + Fu.shape[1]] = Fu

        # form h vector
        for i in range(0, self.N-1):
            if k + i < self.opt.nk - 1:
                u = self.opt.u.value[:, [k + i]]
                indx_h = i*Fu.shape[0]
                
                h_i = np.array([tmax-u[0, 0],tmax+u[0, 0],tmax-u[1, 0],tmax+u[1, 0],tmax-u[2, 0],tmax+u[2, 0]])#tmin * (u.T@u)[0,0]
                self.h[indx_h : indx_h + Fu.shape[0]] = h_i
            
        self.Q.fill(0)
        # form Q matrix
        for i in range(0, self.N-1):
            indx_Q = i*(14 + 3)
            if k + i < self.N-1:
            
                self.Q[indx_Q: indx_Q + 3, indx_Q : indx_Q + 3] = self.F
                self.Q[indx_Q + 3 : indx_Q + 17, indx_Q + 3 : indx_Q + 17] = self.H
            else:
                self.Q[indx_Q: indx_Q + 3, indx_Q : indx_Q + 3] = self.F
                self.Q[indx_Q + 3 : indx_Q + 17, indx_Q + 3 : indx_Q + 17] = self.H2

    def solve_custom(self, k): 
            qp = qp_module.QP(self.Q, self.q, self.G, self.h, self.A, self.b)
            qp.solve()
            print("custom: " + str(qp.x[0:3]))
            return qp.x[0:3]