import numpy as np
import matplotlib.pyplot as plt


def load_data():
    data = {
        "global_s": [],
        "global_v": [],
        "global_e": [],
        "local_s": []
    }

    for relaxation_factor in w_global:
        with open(f'global_error_{relaxation_factor}.txt', 'r') as f:
            error_data = []
            for line in f:
                error_data.append([float(x) for x in line.strip().split()])
            data["global_e"].append(np.transpose(np.array(error_data)))

        with open(f'global_{relaxation_factor}_solution.txt', 'r') as f:
            s_data = []
            for line in f:
                s_data.append([float(x) for x in line.strip().split()])
            data["global_v"].append(np.array(s_data))

        with open(f'global_{relaxation_factor}.txt', 'r') as f:
            s_data = []
            for line in f:
                s_data.append(float(line))
            data["global_s"].append(np.array(s_data))

    for relaxation_factor in w_local:
        with open(f'local_{relaxation_factor}.txt', 'r') as f:
            s_data = []
            for line in f:
                s_data.append(float(line))
            data["local_s"].append(np.array(s_data))

    return data


def plot_relaxation(data, w_global, w_local):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

    for idx, label in enumerate(w_global):
        iterations = list(range(1, len(data["global_s"][idx]) + 1))
        ax1.plot(iterations, data["global_s"][idx], label=f"w={'.'.join(label.split('_'))}")
    ax1.set_xlabel("Iteration number")
    ax1.set_ylabel("S")
    ax1.set_xlim([2, 100000])
    ax1.set_ylim([0, 5000])
    ax1.set_title("Global Relaxation")
    ax1.set_xscale("log")
    ax1.legend()

    for idx, label in enumerate(w_local):
        iterations = list(range(1, len(data["local_s"][idx]) + 1))
        ax2.plot(iterations, data["local_s"][idx], label=f"w={'.'.join(label.split('_'))}")
    ax2.set_xlabel("Iteration number")
    ax2.set_ylabel("S")
    ax2.set_xlim([2, 100000])
    ax2.set_ylim([0, 3000])
    ax2.set_title("Local Relaxation")
    ax2.set_xscale("log")
    ax2.legend()

    plt.savefig("relaxation.png")


def plot_global_results(data):
    fig, axes = plt.subplots(2, 2, figsize=(19, 14))

    plot_grid(axes[0, 0], data["global_v"][0], "Potential: Global Relaxation w=0.6", 150, 100)
    plot_grid(axes[1, 0], data["global_e"][0], "Error: Global Relaxation w=0.6", 150, 100)

    plot_grid(axes[0, 1], data["global_v"][1], "Potential: Global Relaxation w=1.0", 150, 100)
    plot_grid(axes[1, 1], data["global_e"][1], "Error: Global Relaxation w=1.0", 150, 100)

    plt.savefig("global_details.png")


def plot_grid(ax, grid_data, title, x_lim, y_lim):
    im = ax.imshow(grid_data, cmap="plasma", origin="lower")
    ax.set_xlim(0, x_lim)
    ax.set_ylim(0, y_lim)
    ax.set_title(title)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    plt.colorbar(im, ax=ax)


if __name__ == "__main__":
    w_global = ["0_6", "1_0"]
    w_local = ["1_0", "1_4", "1_8", "1_9"]

    data = load_data()

    plot_relaxation(data, w_global, w_local)

    plot_global_results(data)
