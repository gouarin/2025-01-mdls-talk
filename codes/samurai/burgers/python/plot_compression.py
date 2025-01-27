import numpy as np
import h5py
import matplotlib.pyplot as plt
from matplotlib import animation

# color_palette = [
#     "#E63946",
#     "#F1FAEE",
#     "#A8DADC",
#     "#457B9D",
#     "#1D3557"
# ]

color_palette = [
    "#D32F2F",
    "#388E3C",
    "#FBC02D",
    "#7B1FA2",
]


def line_plot(ax, x, y, c, m="o-", alpha=1, linewidth=1.5):
    return ax.plot(x, y, m, color=c, linewidth=linewidth, markersize=1, alpha=alpha)[0]


def read_h5file(filename):
    return h5py.File(filename + ".h5", "r")["mesh"]


def get_mesh_data(h5file):
    points = h5file["points"]
    connectivity = h5file["connectivity"]
    segments = np.zeros((connectivity.shape[0], 2, 2))
    segments[:, :, 0] = points[:][connectivity[:]][:, :, 0]
    centers = 0.5 * (segments[:, 0, 0] + segments[:, 1, 0])
    index = np.argsort(centers)
    return centers, index


def get_field_data(h5file, name):
    return h5file["fields"][name][:]


def ax_properties(ax, fontsize=4):
    for a in ax.flatten():
        for label in a.get_xticklabels() + a.get_yticklabels():
            label.set_fontsize(fontsize)
        # a.set_facecolor('#F1FAEE30')
        # a.grid(True, color='#A8DADC', axis='y')
        a.grid(True, color="#9E9E9E", axis="y")
        a.yaxis.set_major_locator(
            plt.MaxNLocator(6)
        )  # Set number of ticks in y direction

    if ax.ndim == 2:
        for a in ax:
            a[-1].yaxis.tick_right()  # Put the y-axis ticks on the right
            for label in a[-1].get_xticklabels() + a[-1].get_yticklabels():
                label.set_fontsize(fontsize)
    else:
        ax[-1].yaxis.tick_right()  # Put the y-axis ticks on the right
        for label in ax[-1].get_xticklabels() + ax[-1].get_yticklabels():
            label.set_fontsize(fontsize)


fig = plt.figure()
ax = fig.subplots(4, 3)

ax_properties(ax)

test_cases = ["build/exp", "build/abs", "build/sqrt", "build/tanh"]
test_names = [
    "$exp(-50x^2)$",
    "$1 - |2x|$ if $-0.5 < x < 0.5$, $0$ elsewhere",
    r"$1 - \sqrt{\left| sin \left( \frac{\pi}{2} x \right) \right|}$",
    "$tanh(50 |x|) - 1$",
]

linewidth = 0.75

ind = 0
for c in test_cases:
    data = read_h5file(c)
    centers, index = get_mesh_data(data)
    u = get_field_data(data, "u")
    level = get_field_data(data, "level")
    error = get_field_data(data, "error")
    ax[ind, 1].set_ylim([1, 13])

    ax[ind, 0].text(
        0.05,
        0.92,
        f"{test_names[ind]}",
        transform=ax[ind, 0].transAxes,
        fontsize=4,
        verticalalignment="top",
        usetex=True,
    )

    line_plot(
        ax[ind, 0],
        centers[index],
        u[index],
        color_palette[1],
        m="",
        linewidth=linewidth,
    )
    line_plot(
        ax[ind, 1],
        centers[index],
        level[index],
        color_palette[0],
        m="",
        linewidth=linewidth,
    )
    line_plot(
        ax[ind, 2],
        centers[index],
        error[index],
        color_palette[3],
        m="",
        linewidth=linewidth,
    )
    ind += 1
plt.savefig("compression.png", dpi=300)


linewidth = 1.5
ind = 0
for c in test_cases:
    fig = plt.figure(figsize=(12, 4))
    ax = fig.subplots(1, 3)
    ax_properties(ax, 8)

    data = read_h5file(c)
    centers, index = get_mesh_data(data)
    u = get_field_data(data, "u")
    level = get_field_data(data, "level")
    error = get_field_data(data, "error")
    ax[1].set_ylim([1, 13])

    ax[0].set_title("solution")
    line_plot(
        ax[0],
        centers[index],
        u[index],
        color_palette[1],
        m="",
        linewidth=linewidth,
    )

    ax[1].set_title("level")
    line_plot(
        ax[1],
        centers[index],
        level[index],
        color_palette[0],
        m="",
        linewidth=linewidth,
    )

    ax[2].set_title("error")
    line_plot(
        ax[2],
        centers[index],
        error[index],
        color_palette[3],
        m="",
        linewidth=linewidth,
    )
    ind += 1
    plt.savefig(f"compression_{c.split('/')[-1]}.png", dpi=300)
