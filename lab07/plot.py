import numpy as np
import matplotlib.pyplot as plt


def load_grid(filename):
    with open(filename, "r") as file:
        return [float(value) for value in file.readlines()]


def load_data(filename):
    data = []
    with open(filename, "r") as file:
        for line in file.readlines():
            data.append([float(value) for value in line.split()])
    return np.array(data)


# Wczytywanie siatek
x_grid = load_grid("data/x.txt")
y_grid = load_grid("data/y.txt")

# Lista wartości Q i odpowiadających im plików
Q_values = [-1000, -4000, 4000]
files = {
    -1000: {
        "psi": "data/psi_-1000.txt",
        "zeta": "data/zeta_-1000.txt",
        "u": "data/u_-1000.txt",
        "v": "data/v_-1000.txt",
        "clim_psi": [-55, -50],
        "clim_zeta": [-200, 300]
    },
    -4000: {
        "psi": "data/psi_-4000.txt",
        "zeta": "data/zeta_-4000.txt",
        "u": "data/u_-4000.txt",
        "v": "data/v_-4000.txt",
        "clim_psi": [-218, -202],
        "clim_zeta": [-700, 1100]
    },
    4000: {
        "psi": "data/psi_4000.txt",
        "zeta": "data/zeta_4000.txt",
        "u": "data/u_4000.txt",
        "v": "data/v_4000.txt",
        "clim_psi": [202, 218],
        "clim_zeta": [-200, 350]
    },

}

# Generowanie wykresów dla każdej wartości Q
for Q in Q_values:
    psi = np.transpose(load_data(files[Q]["psi"]))
    zeta = np.transpose(load_data(files[Q]["zeta"]))
    u = np.transpose(load_data(files[Q]["u"]))
    v = np.transpose(load_data(files[Q]["v"]))

    # Rysowanie wykresów
    fig, axs = plt.subplots(2, 2, figsize=(19, 14))

    # Wykres 1: Kontury psi
    cs1 = axs[0, 0].contour(psi, cmap='inferno', levels=500)
    cs1.set_clim(files[Q]["clim_psi"][0], files[Q]["clim_psi"][1])
    axs[0, 0].set_xlim(0, 200)
    axs[0, 0].set_ylim(0, 90)
    axs[0, 0].set_xlabel('x')
    axs[0, 0].set_ylabel('y')
    axs[0, 0].title.set_text(f"Q={Q}, ψ(x,y)")
    axs[0, 0].set_aspect('equal', adjustable='box')

    # Wykres 2: Kontury zeta
    cs2 = axs[0, 1].contour(zeta, cmap='inferno', levels=60)
    cs2.set_clim(files[Q]["clim_zeta"][0], files[Q]["clim_zeta"][1])
    axs[0, 1].set_xlim(0, 200)
    axs[0, 1].set_ylim(0, 90)
    axs[0, 1].set_xlabel('x')
    axs[0, 1].set_ylabel('y')
    axs[0, 1].title.set_text(f"Q={Q}, ζ(x,y)")
    axs[0, 1].set_aspect('equal', adjustable='box')

    # Wykres 3: Prędkość u
    cs3 = axs[1, 0].imshow(u, cmap="jet", origin="lower", extent=[min(x_grid), max(x_grid), min(y_grid), max(y_grid)])
    fig.colorbar(cs3, ax=axs[1, 0], fraction=0.022, pad=0.01)
    axs[1, 0].set_title(f"Q={Q}, u(x,y)")
    axs[1, 0].set_xlabel("x")
    axs[1, 0].set_ylabel("y")

    # Wykres 4: Prędkość v
    cs4 = axs[1, 1].imshow(v, cmap="jet", origin="lower", extent=[min(x_grid), max(x_grid), min(y_grid), max(y_grid)])
    fig.colorbar(cs4, ax=axs[1, 1], fraction=0.022, pad=0.01)
    axs[1, 1].set_title(f"Q={Q}, v(x,y)")
    axs[1, 1].set_xlabel("x")
    axs[1, 1].set_ylabel("y")

    plt.tight_layout()
    plt.savefig(f"contour_plots_Q{Q}.png", bbox_inches='tight', transparent=False)
    plt.close()
