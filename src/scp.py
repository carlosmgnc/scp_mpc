import numpy as np
import cvxpy as cvx
import sympy as sy
from sympy import *
from cvxpygen import cpg
from scp_socp_solver.cpg_solver import cpg_solve 
import time

class optProblem:
    def __init__(self):

        # vehicle properties
        self.mw = 2.0
        self.md = 1.0
        self.Tmax = 5
        self.Tmin = 1
        self.rt = np.array([[-0.01], [0], [0]])
        self.Jbvec = np.array([[0.01], [0.01], [0.01]])
        self.Jb = np.diag(self.Jbvec.flatten())
        self.Jbinv = np.diag((1 / self.Jbvec).flatten())
        self.alpha = 0.07

        #discrete time grid
        self.nk = 35
        self.K = np.arange(0, self.nk)
        self.dt = 1 / (self.nk - 1)
        self.tau = np.linspace(0, 1, self.nk)

        # initial trajectory guess
        ri = np.array([[4], [2], [0]])
        vi = np.array([[-1], [-1], [1]])
        self.vf = np.array([[0], [0], [0]])
        self.qi = np.array([[1],[0],[0],[0]])
        self.g = np.array([[-1], [0], [0]])

        self.tfguess = 5
        alpha1_list = (1 - self.tau) / (1)
        alpha2_list = 1 - alpha1_list

        self.mk = alpha1_list * self.mw + alpha2_list * self.md
        self.rk = alpha1_list * ri
        self.vk = alpha1_list * vi + alpha2_list * self.vf
        
        self.qk = np.tile(self.qi, (1, self.nk))
        self.wk = np.zeros((3, self.nk))

        # optimization problem
        self.x = cvx.Variable((14, self.nk))
        self.u = cvx.Variable((3, self.nk))
        self.sigma = cvx.Variable(nonneg=True)
        self.nu = cvx.Variable((14, self.nk - 1))
        self.delta = cvx.Variable((self.nk, 1), nonneg=True)
        self.delta_sigma = cvx.Variable(nonneg=True)

        self.xk = cvx.Parameter((14, self.nk))
        self.uk = cvx.Parameter((3, self.nk))
        self.sigmak = cvx.Parameter(nonneg=True)
        self.A = cvx.Parameter((14, 14 * (self.nk - 1)))
        self.B = cvx.Parameter((14, 3 * (self.nk - 1)))
        self.C = cvx.Parameter((14, 3 * (self.nk - 1)))
        self.E = cvx.Parameter((14, self.nk - 1))
        self.z = cvx.Parameter((14, self.nk - 1))
        self.w_nu = cvx.Parameter(nonneg=True)
        self.w_delta = cvx.Parameter(nonneg=True)
        self.w_sigma = cvx.Parameter(nonneg=True)

        self.cos_theta_max = cvx.Parameter(nonneg=True)
        self.cos_delta_max = cvx.Parameter(nonneg=True)
        self.w_max = cvx.Parameter(nonneg=True)
        self.tan_gamma_gs = cvx.Parameter(nonneg=True)
        
        self.cos_theta_max.value = np.cos(np.deg2rad(30))
        self.cos_delta_max.value = np.cos(np.deg2rad(20))
        self.w_max.value = np.deg2rad(90)
        self.tan_gamma_gs.value = np.tan(np.deg2rad(10))
        self.w_nu.value = 10000
        self.w_delta.value = 1
        self.w_sigma.value = 0.1

        # initialize problem parameters

        self.xk.value = np.vstack([self.mk, self.rk, self.vk, self.qk, self.wk])
        self.uk.value = -self.mk * self.g
        self.sigmak.value = self.tfguess

        self.A.value = np.zeros((14, 14 * (self.nk - 1)))
        self.B.value = np.zeros((14, 3 * (self.nk - 1)))
        self.C.value = np.zeros((14, 3 * (self.nk - 1)))
        self.E.value = np.zeros((14, self.nk - 1))
        self.z.value = np.zeros((14, self.nk - 1))
        self.B_rcs = np.zeros((14, self.nk - 1))
        self.B_rcs_ct = np.zeros((14, 1))
        self.B_rcs_ct[11, 0] = 1*self.Jbinv[0,0]

        print("b_rcs: " + str(self.B_rcs_ct))

        self.stm_inv = np.zeros((14, 14))
        self.A_f, self.B_f, self.E_f, self.z_f = self.def_jacobian_funcs()

        self.constraints = []
        self.cost = 0

    # Direction Cosine Matrix Function
    def DCM(self, q): 
        return sy.Matrix(
            [
                [
                    1 - 2 * (q[2] ** 2 + q[3] ** 2),
                    2 * (q[1] * q[2] + q[0] * q[3]),
                    2 * (q[1] * q[3] - q[0] * q[2]),
                ],
                [
                    2 * (q[1] * q[2] - q[0] * q[3]),
                    1 - 2 * (q[1] ** 2 + q[3] ** 2),
                    2 * (q[2] * q[3] + q[0] * q[1]),
                ],
                [
                    2 * (q[1] * q[3] + q[0] * q[2]),
                    2 * (q[2] * q[3] - q[0] * q[1]),
                    1 - 2 * (q[1] ** 2 + q[2] ** 2),
                ],
            ]
        )
    
    # skew symmetric quaternion matrix
    def omega(self, w):
        return sy.Matrix(
        [
            [0, -w[0], -w[1], -w[2]],
            [w[0], 0, w[2], -w[1]],
            [w[1], -w[2], 0, w[0]],
            [w[2], w[1], -w[0], 0],
        ]
    )

    # skew symmetric cross product matrix function
    def cr(self, v):
        return sy.Matrix([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    
    # returns linear interpolation of a vector between two time steps 
    def get_alphas(self, t):
        k = int(np.floor(t / self.dt))
        tk = self.tau[k]

        if k + 1 <= self.nk-1:
            tk1 = self.tau[k + 1]
        else: 
            tk1 = 1 + self.dt

        alphak = (tk1 - t) / (tk1 - tk)
        betak = 1 - alphak

        return alphak, betak, tk, tk1, k
    
    # returns linearized system matrices as a function
    def def_jacobian_funcs(self):
        f_symb = sy.zeros(14, 1)
        x = sy.Matrix(sy.symbols("m r0 r1 r2 v0 v1 v2 q0 q1 q2 q3 w0 w1 w2", real=True))
        u = sy.Matrix(sy.symbols("u0 u1 u2", real=True))

        g = sy.Matrix(self.g)
        rt = sy.Matrix(self.rt)
        Jb = sy.Matrix(self.Jb)
        Jbinv = sy.Matrix(self.Jbinv)

        f_symb[0, 0] = -self.alpha * u.norm()
        f_symb[1:4, 0] = x[4:7, 0]
        f_symb[4:7, 0] = (1/x[0, 0]) * self.DCM(x[7:11, 0]).T * u + g
        f_symb[7:11, 0] = (1/2) * self.omega(x[11:14, 0]) * x[7:11, 0]
        f_symb[11:14, 0] = Jbinv * (self.cr(rt) * u - self.cr(x[11:14, 0]) * Jb * x[11:14, 0])

        f_symb = f_symb
        A_symb = f_symb.jacobian(x) * self.sigmak.value
        B_symb = f_symb.jacobian(u) * self.sigmak.value
        z_symb = - A_symb * x - B_symb * u

        A_f = sy.lambdify((x, u), A_symb, 'numpy')
        B_f = sy.lambdify((x, u), B_symb, 'numpy')
        E_f = sy.lambdify((x, u), f_symb, 'numpy')
        z_f = sy.lambdify((x, u), z_symb, 'numpy')

        return A_f, B_f, E_f, z_f

    # derivative function for rk4
    def P_dot(self, t, X):
        phi = X[14 : 14 + 14 * 14].reshape((14, 14))
        phi_inv = np.linalg.inv(phi)

        alphak, betak, tk, tk1, k = self.get_alphas(t)

        if k + 1 <= self.nk-1:
            u = (alphak * self.uk.value[:, k] + betak * self.uk.value[:, k+1]).reshape((3, 1))
        else: 
            u = (self.uk.value[:, k]).reshape((3, 1))

        X_flat = X[:14].flatten()
        u_flat = u.flatten()

        P1 = self.E_f(X_flat, u_flat) * self.sigmak.value
        P2 = (self.A_f(X_flat, u_flat) @ phi).reshape((14 * 14, 1))
        P3 = (phi_inv @ self.B_f(X_flat, u_flat) * alphak).reshape((14 * 3, 1))
        P4 = (phi_inv @ self.B_f(X_flat, u_flat) * betak).reshape((14 * 3, 1))
        P5 = phi_inv @ self.E_f(X_flat, u_flat)
        P6 = phi_inv @ self.z_f(X_flat, u_flat)

        P_RCS = self.sigmak.value * phi_inv @ self.B_rcs_ct

        return np.vstack([P1, P2, P3, P4, P5, P6, P_RCS])

    # rk4 single step function
    def rk41(self, func, tk, xk, dt):
        k1 = func(tk, xk)
        k2 = func(tk + dt / 2, xk + (dt / 2) * k1)
        k3 = func(tk + dt / 2, xk + (dt / 2) * k2)
        k4 = func(tk + dt, xk + dt * k3)
        output = xk + (dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
        return output

    # discretization using multiple shooting
    def discretize(self):
        nsub = 2
        dt_sub = self.dt / (nsub + 1)

        #indeces for flattened state
        n = 14
        m = 3
        nsq = n * n
        nm = n * m

        xnd = n
        stmnd = xnd + nsq
        Bnd = stmnd + nm
        Cnd = Bnd + nm
        End = Cnd + n
        znd = End + n
        rcsnd = znd + n

        P_temp = np.empty((rcsnd, 1))

        for i in range(0, self.nk - 1):
            P_temp[:xnd, 0] = self.xk.value[:, i]
            P_temp[xnd: stmnd, 0] = np.eye(n).reshape(nsq)
            P_temp[stmnd : Bnd, 0] = np.zeros(nm)
            P_temp[Bnd: Cnd, 0] = np.zeros(nm)
            P_temp[Cnd: End, 0] = np.zeros(n)
            P_temp[End : znd, 0] = np.zeros(n)
            P_temp[znd : rcsnd, 0] = np.zeros(n)

            for j in range(0, nsub+1):
                sub_time = i * self.dt + j * dt_sub
                P_temp = self.rk41(self.P_dot, sub_time, P_temp, dt_sub)

            stm = P_temp[xnd : stmnd].reshape((n, n))
            indx1 = i * n
            indx2 = i * m

            self.A.value[:, indx1:indx1 + 14] = stm
            self.B.value[:, indx2:indx2 + 3] = stm @ P_temp[stmnd : Bnd].reshape((n, m))
            self.C.value[:, indx2:indx2 + 3] = stm @ P_temp[Bnd : Cnd].reshape((n, m))
            self.E.value[:, i] = stm @ P_temp[Cnd : End, 0]
            self.z.value[:, i] = stm @ P_temp[End : znd, 0]
            self.B_rcs[:, i] = stm @ P_temp[znd: rcsnd, 0]

    def def_cvx_problem(self):
        nu_cost = 0
        self.constraints += [self.md <= self.x[0, :]]

        #boundary conditions
        self.constraints += [self.x[0, 0] == self.xk[0, 0]]
        self.constraints += [self.x[1:4, 0] == self.xk[1:4, 0]]
        self.constraints += [self.x[4:7, 0] == self.xk[4:7, 0]]
        self.constraints += [self.x[11:14, 0] == self.xk[11:14, 0]]

        self.constraints += [self.x[1:, -1] == self.xk[1:, -1]]
        self.constraints += [self.u[1:3, -1] == np.zeros(2)]

        e2 = np.array([[0], [1], [0]])
        e3 = np.array([[0], [0], [1]])
        H23 = np.hstack([e2, e3])

        for k in range(self.nk):
            # difference in state, control and time from last iterate
            dx = self.x[:,[k]]-self.xk[:, [k]]
            du = self.u[:, [k]]-self.uk[:, [k]]
            dsigma = self.sigma - self.sigmak

            #state and control constraints
            H = np.array([[0, 0, 1, 0], [0, 0, 0, 1]])
            self.constraints += [self.cos_theta_max <= 1 - 2 * cvx.square(cvx.norm(H @ self.x[7:11, [k]]))]
            self.constraints += [cvx.norm(self.x[11:14, [k]]) <= self.w_max]
            self.constraints += [cvx.norm(self.u[:, [k]]) <= self.Tmax]
            self.constraints += [self.Tmin*cvx.norm(self.uk[:, [k]]) <= cvx.transpose(self.uk[:, [k]]) @ self.u[:, [k]]]
            self.constraints += [self.cos_delta_max * cvx.norm(self.u[:, [k]]) <= self.u[0, [k]]]

            self.constraints += [self.tan_gamma_gs * cvx.norm(H23.T @ self.x[1:4, [k]]) <= self.x[1, [k]]]
            self.constraints += [cvx.square(cvx.norm(dx)) + cvx.square(cvx.norm(du)) <= self.delta[k,0] ]
            self.constraints += [cvx.norm(dsigma, 1) <= self.delta_sigma]
        
        #dynamics constrains
        for k in range(0, self.nk - 1):
            indx1 = k*14
            indx2 = k*3
            self.constraints += [
                self.x[:, [k + 1]]
                == self.A[:, indx1:indx1 + 14] @ self.x[:, [k]]
                + self.B[:, indx2:indx2 + 3] @ self.u[:, [k]]
                + self.C[:, indx2:indx2 + 3] @ self.u[:, [k + 1]]
                + self.E[:, [k]] * self.sigma
                + self.z[:, [k]]
                + self.nu[:, [k]]
            ]

            nu_cost += cvx.norm(self.nu[:, [k]])
        
        self.cost = (
            self.sigma
            + self.w_nu * nu_cost
            + self.w_delta * cvx.norm(self.delta)
            + self.w_sigma * cvx.norm(self.delta_sigma)
        )

        objective = cvx.Minimize(self.cost)
        self.prob = cvx.Problem(objective, self.constraints)

        #cpg.generate_code(self.prob, code_dir='scp_socp_solver', solver='QOCO', prefix="scp_socp_")
        self.prob.register_solve("scp_socp", cpg_solve) 
    
    def solve_cvx_problem(self):
        print("solving")
        print("----------------------------------")
        #self.prob.solve(solver="QOCO", ignore_dpp=True)
        self.prob.solve(method="scp_socp")
        print("solver status : " + self.prob.status)
        print("solve time    : " + f"{self.prob.solver_stats.solve_time * 1000:.2f}" + " ms")
        print("cost          : " + f"{self.prob.objective.value:.3f}")
        print("sigma         : " + f"{self.sigma.value:.3f}")
        print("----------------------------------")

        self.xk.value = self.x.value
        self.uk.value = self.u.value
        self.sigmak.value = self.sigma.value

        return self.x.value, self.u.value, self.sigma.value, self.cost.value, self.nu.value 