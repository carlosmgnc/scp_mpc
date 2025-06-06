import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.animation import FuncAnimation

class Plots():

    def __init__(self, anim, skip):
        self.fig_traj_plot = plt.figure(7, figsize=(8, 8))
        self.ax = plt.axes(projection="3d")
        self.anim = anim
        self.skip = skip
        
    def DCM_output(self, q): 
        return np.array(
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

    def plot(self, solver, sim):

        x = sim.trajectory
        opt = solver.opt
        u = solver.u
        nu = solver.nu

        plt.figure(1)
        plt.title("pos vs time")
        labels = []

        for i in range(3):
            plt.plot(opt.tau, solver.x[1 + i, :], label="", linestyle="--")

        for i in range(3):
            plt.plot(opt.tau, x[1 + i, :], label="")
        plt.legend(["x", "y", "z"])
        plt.xlabel("time")
        plt.ylabel("position")


        plt.figure(2)
        plt.title("mass vs time")

        plt.plot(opt.tau, x[0, :], label="")
        plt.xlabel("time")
        plt.ylabel("mass")

        plt.figure(3)
        plt.title("control input vs time")
        labels = []

        for i in range(3):
            plt.plot(opt.tau, u[i, :], label="", linestyle="--")

        for i in range(3):
            plt.plot(opt.tau, u[i, :] + sim.u_mpc_list[i, :], label="")

        plt.legend(["ux", "uy", "uz"])
        plt.xlabel("time")
        plt.ylabel("thrust vector")

        plt.figure(4)
        plt.title("$\| T(t) \|$  vs Time")
        plt.plot(opt.tau[:], np.linalg.norm(u[:, :] + sim.u_mpc_list[:, :], axis=0), color="lightgrey")
        plt.plot(opt.tau[:], np.linalg.norm(u[:, :], axis=0), linestyle="--", color="blue")
        plt.plot([opt.tau[0], opt.tau[-1]], [opt.Tmax, opt.Tmax], linestyle='--', color='black', linewidth=1)
        plt.plot([opt.tau[0], opt.tau[-1]], [opt.Tmin, opt.Tmin], linestyle='--', color='black', linewidth=1)
        plt.xlabel("time")
        plt.ylabel("$\| T(t) \|$")

        plt.figure(5)
        plt.title("virtual control")
        labels = []

        for i in range(nu.shape[0]):
            plt.plot(opt.tau[1:], nu[i, :])
            #print(np.linalg.norm(max(nu[i, :])))

            plt.xlabel("time")
            plt.ylabel("virtual control")
        plt.legend()

        # convergence plot
        plt.figure(6)
        plt.title("final time vs iteration number")
        plt.plot(range(solver.converged_iter), solver.sigma_list[:solver.converged_iter])
        plt.xlabel("iteration")
        plt.ylabel("final time")

        # 3d trajectory plot
        # fig_traj_plot = plt.figure(7, figsize=(8, 8))
        # fig_traj_plot.subplots_adjust(left=0.1, right=0.9, top=0.95, bottom=0.05)
        self.fig_traj_plot.tight_layout()
        # ax = plt.axes(projection="3d")
        self.ax.view_init(elev=25, azim=161)

        #ax.view_init(elev=0, azim=90)
        self.ax.plot3D(x[2, :], x[3, :], x[1, :], color="darkgrey")

        # fix aspect ratio of 3d plot
        x_lim = self.ax.get_xlim3d()
        y_lim = self.ax.get_ylim3d()
        z_lim = self.ax.get_zlim3d()

        max_lim = max(
            abs(x_lim[1] - x_lim[0]), abs(y_lim[1] - y_lim[0]), abs(z_lim[1] - z_lim[0])
        )
        x_mid = sum(x_lim) * 0.5
        y_mid = sum(y_lim) * 0.5

        rt_I = np.zeros((3, opt.nk))
        for i in range(opt.nk):
            rt_I[:, [i]] = 15 * self.DCM_output(x[7:11, i]).T @ opt.rt

        thrust_vecs = np.empty((3, opt.nk))
        qlen = 0.03 * max_lim

        for i in range(rt_I.shape[1]):
            thrust_vecs[:, [i]] = self.DCM_output(x[7:11, i]).T @ (u[:, [i]] + sim.u_mpc_list[:, [i]])

        q = qlen * thrust_vecs

        base_x = x[1, :] - q[0, :]
        base_y = x[2, :] - q[1, :]
        base_z = x[3, :] - q[2, :]

        self.ax.quiver(
            base_y[::self.skip],
            base_z[::self.skip],
            base_x[::self.skip],
            q[1, ::self.skip],
            q[2, ::self.skip],
            q[0, ::self.skip],
            normalize=False,
            arrow_length_ratio=0.1,
            # color=(1, 60/255, 0),
            linewidth=1,
        )

        self.ax.quiver(
            base_y[-1],
            base_z[-1],
            base_x[-1],
            q[1, -1],
            q[2, -1],
            q[0, -1],
            normalize=False,
            arrow_length_ratio=0.1,
            # color=(1, 60/255, 0),
            linewidth=1,
        )

        base_x_2 = x[1, :]
        base_y_2 = x[2, :]
        base_z_2 = x[3, :]

        self.ax.quiver(
            base_y_2[::self.skip],
            base_z_2[::self.skip],
            base_x_2[::self.skip],
            -2 * rt_I[1, ::self.skip],
            -2 * rt_I[2, ::self.skip],
            -2 * rt_I[0, ::self.skip],
            normalize=False,
            arrow_length_ratio=0,
            color=(0.1, 0.1, 0.1),
            linewidth=2.0,
        )

        self.ax.quiver(
            base_y_2[-1],
            base_z_2[-1],
            base_x_2[-1],
            -2 * rt_I[1, -1],
            -2 * rt_I[2, -1],
            -2 * rt_I[0, -1],
            normalize=False,
            arrow_length_ratio=0,
            color=(0.1, 0.1, 0.1),
            linewidth=2.0,
        )

        self.ax.set_xlim3d([x_mid - max_lim * 0.5, x_mid + max_lim * 0.5])
        self.ax.set_ylim3d([y_mid - max_lim * 0.5, y_mid + max_lim * 0.5])
        self.ax.set_zlim3d([0, max_lim])
        self.ax.plot(self.ax.get_xlim(), (0, 0), (0, 0), color="black", linestyle="--", linewidth=1)
        self.ax.plot((0, 0), self.ax.get_ylim(), (0, 0), color="black", linestyle="--", linewidth=1)
        # self.ax.plot(solver.trajectory_list[-1, 2, :], solver.trajectory_list[-1, 3, :], solver.trajectory_list[-1, 1, :], linestyle="--", linewidth=0.5, color="black")
        self.ax.plot(opt.x.value[2, :], opt.x.value[3, :], opt.x.value[1, :], linestyle="--", linewidth=2, color="blue")


        def shared_traj_plot_properties(ax):
            ax.set_title("converged trajectory")
            ax.scatter(0, 0, 0, color="green", s=10)
            ax.set_xlabel("y")
            ax.set_ylabel("z")
            ax.set_zlabel("x")


        shared_traj_plot_properties(self.ax)

        ############################# animation #############################
        matplotlib.rcParams['axes3d.mouserotationstyle'] = 'azel'

        if self.anim == True:
            fig_anim = plt.figure(8, figsize=(8, 8))
            fig_anim.tight_layout()
            ax_anim = plt.axes(projection="3d")
            # ax_anim.view_init(elev=25, azim=161)
            ax_anim.view_init(elev=11, azim=-135)
            ax_anim.plot3D(x[2, :], x[3, :], x[1, :], linestyle="--", linewidth=0.5, color="black")
            #solver.converged_iter
            # for i in range(solver.converged_iter):
            #     #print(trajectory_list[i, 2, :])
            ax_anim.plot3D(solver.trajectory_list[-1, 2, :], solver.trajectory_list[-1, 3, :], solver.trajectory_list[-1, 1, :], linestyle="--", linewidth=0.5, color="blue")

            shared_traj_plot_properties(ax_anim)
            ax_anim.set_xlim(self.ax.get_xlim())
            ax_anim.set_ylim(self.ax.get_ylim())
            ax_anim.set_zlim(self.ax.get_zlim())

            quiver = ax_anim.quiver(
                base_y[0],
                base_z[0],
                base_x[0],
                q[1, 0],
                q[2, 0],
                q[0, 0],
                normalize=False,
                arrow_length_ratio=0.1,
                color=(1, 60/255, 0),
                linewidth=1,
            )
            quiver2 = ax_anim.quiver(
                base_y_2[0],
                base_z_2[0],
                base_x_2[0],
                -2 * rt_I[1, 0],
                -2 * rt_I[2, 0],
                -2 * rt_I[0, 0],
                normalize=False,
                arrow_length_ratio=0,
                color=(0.1, 0.1, 0.1),
                linewidth=2.0,
            )


            def update(frame):
                quiver.set_segments(
                    [
                        [
                            [base_y[frame], base_z[frame], base_x[frame]],
                            [
                                x[2, frame],
                                x[3, frame],
                                x[1, frame],
                            ],
                        ]
                    ]
                )

                quiver2.set_segments(
                    [
                        [
                            [base_y_2[frame], base_z_2[frame], base_x_2[frame]],
                            [
                                base_y_2[frame] - 2 * rt_I[1, frame],
                                base_z_2[frame] - 2 * rt_I[2, frame],
                                base_x_2[frame] - 2 * rt_I[0, frame],
                            ],
                        ]
                    ]
                )

                return quiver, quiver2


            anim_int = 100
            animation = FuncAnimation(fig_anim, update, frames=opt.nk, interval=anim_int)

            fig_names = ["position", "mass", "control", "throttle", "virtual_control", "tof_iteration", "trajectory", "animation"]

            # for i in range(1, 8):
            #     plt.figure(i).savefig("../images/" + fig_names[i - 1] + ".png", dpi=300)

            animation.save("../images/animation.gif", writer="pillow", fps=1000 / anim_int)

            plt.show(block=False)
            plt.pause(1)
            input()
            plt.close()
        
        # plt.figure()
        # plt.plot(opt.tau, sim.u_rcs_list)
        # plt.title("rcs control")

        # plt.figure()
        # plt.title("quat vs time")
        # labels = []

        # for i in range(4):
        #     plt.plot(opt.tau, solver.x[7 + i, :], label="", linestyle="--")

        # for i in range(4):
        #     plt.plot(opt.tau, x[7 + i, :], label="")
        # plt.legend(["w", "qx", "qy", "qz", "w", "qx", "qy", "qz"])
        # plt.xlabel("time")
        # plt.ylabel("position")
        
