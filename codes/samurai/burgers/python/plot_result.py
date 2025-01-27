import numpy as np
import h5py
import matplotlib.pyplot as plt
from matplotlib import animation

ref = "./build/finest/burgers_1D_ite"


def line_plot(ax, x, y, c, m="o-", alpha=1, linewidth=1.5):
    return ax.plot(x, y, m, color=c, linewidth=linewidth, markersize=2, alpha=alpha)[0]


def read_mesh(filename, ite=None):
    return h5py.File(filename + ".h5", "r")["mesh"]


def get_data(filename):
    mesh = read_mesh(filename)
    points = mesh["points"]
    connectivity = mesh["connectivity"]
    segments = np.zeros((connectivity.shape[0], 2, 2))
    segments[:, :, 0] = points[:][connectivity[:]][:, :, 0]

    u = mesh["fields"]["u"][:]
    if "level" in mesh["fields"]:
        level = mesh["fields"]["level"][:]
    else:
        level = None
    centers = 0.5 * (segments[:, 0, 0] + segments[:, 1, 0])
    index = np.argsort(centers)
    return u, level, centers, index


def make_fig(filename, adapted, adapted_r, lim):
    fig = plt.figure()
    if adapted_r is None:
        ax = fig.subplots(2, 1)
    else:
        ax = fig.subplots(3, 1)

    for a in ax:
        for label in (
            [a.title, a.xaxis.label, a.yaxis.label, a.yaxis.get_offset_text()]
            + a.get_xticklabels()
            + a.get_yticklabels()
        ):
            label.set_fontsize(4)
        a.set_facecolor("#F1FAEE30")
        a.grid(True, color="#A8DADC")

    ax[0].set_ylim([-0.1, 1.1])
    ax[0].text(
        0.01,
        0.95,
        "solution",
        transform=ax[0].transAxes,
        fontsize=10,
        verticalalignment="top",
    )
    ax[1].set_ylim([3.5, 12.5])
    ax[1].text(
        0.01,
        0.95,
        "level",
        transform=ax[1].transAxes,
        fontsize=10,
        verticalalignment="top",
    )
    if adapted_r is not None:
        ax[2].set_ylim([-0.01 * lim, 1.1 * lim])
        ax[2].text(
            0.01,
            0.95,
            "error at the finest level",
            transform=ax[2].transAxes,
            fontsize=10,
            verticalalignment="top",
        )

    start, end = 0, 49

    u, level, centers, index = get_data(f"{adapted}_0")
    line_a = line_plot(ax[0], centers[index], u[index], "#457B9D")
    line_l = line_plot(ax[1], centers[index], level[index], "#E63946")

    u, _, centers, index = get_data(f"{ref}_0")
    line_r = line_plot(
        ax[0], centers[index], u[index], "black", m="", linewidth=1, alpha=0.5
    )

    if adapted_r is not None:
        ur, level, centers, index = get_data(f"{adapted_r}_0")
        # line_ar = line_plot(ax[2], centers[index], np.abs(u[index]- ur[index]), '#1D3557', m='', linewidth=1)
        line_ar = line_plot(
            ax[2],
            centers[index],
            np.abs(u[index] - ur[index]),
            "#1D3557",
            m="",
            linewidth=1,
        )

    def animate(i):
        u, level, centers, index = get_data(f"{adapted}_{i}")
        line_a.set_data(centers[index], u[index])
        line_l.set_data(centers[index], level[index])

        u, _, centers, index = get_data(f"{ref}_{i}")
        line_r.set_data(centers[index], u[index])

        if adapted_r is not None:
            ur, level, centers, index = get_data(f"{adapted_r}_{i}")
            error = np.abs(u[index] - ur[index])
            line_ar.set_data(centers[index], error)
            # ax[-1].set_ylim(error.min() * 1.1, error.max() * 1.1)

    ani = animation.FuncAnimation(fig, animate, frames=end - start, repeat=True)

    writermp4 = animation.FFMpegWriter(fps=1)
    ani.save(f"{filename}.mp4", dpi=300)


min_level = 1
max_level = 12
reg = 0
eps = [0.01, 0.001, 0.0001]

for e in eps:
    name = (
        f"MRA_upwind_without_portion-min_{min_level}-max_{max_level}-eps_{e}-reg_{reg}"
    )
    path = f"./build/{name}/"
    adapted = f"{path}/burgers_1D_ite"
    adapted_r = f"{path}/burgers_1D_recons_ite"
    make_fig(name, adapted, adapted_r, lim=10 * e)

    name = f"MRA_upwind_with_portion-min_{min_level}-max_{max_level}-eps_{e}-reg_{reg}"
    path = f"./build/{name}/"
    adapted = f"{path}/burgers_1D_ite"
    adapted_r = f"{path}/burgers_1D_recons_ite"
    make_fig(name, adapted, adapted_r, lim=e)

# adapted = "./build/MRA_LW_without_portion/burgers_1D_ite"
# adapted_r = "./build/MRA_LW_without_portion/burgers_1D_recons_ite"
# make_fig("LW_without_portion", adapted, adapted_r, lim=1.1)

# adapted_p = "./build/MRA_LW_with_portion/burgers_1D_ite"
# adapted_p_r = "./build/MRA_LW_with_portion/burgers_1D_recons_ite"
# make_fig("LW_with_portion", adapted_p, adapted_p_r, lim=0.04)


# amr_cases = [
#     "AMR_derivative_eps_1e-2_upwind_without_portion",
#     "AMR_derivative_eps_1e-2_upwind_with_portion",
#     "AMR_second_derivative_eps_10_upwind_without_portion",
#     "AMR_second_derivative_eps_10_upwind_with_portion",
#     "AMR_second_derivative_eps_100_upwind_without_portion",
#     "AMR_second_derivative_eps_100_upwind_with_portion",
# ]

# for c in amr_cases:
#     adapted_p = f"./build/{c}/burgers_1D_ite"
#     make_fig(c, adapted_p, None, lim=0.04)
