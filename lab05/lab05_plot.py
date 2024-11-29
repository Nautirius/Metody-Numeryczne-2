import numpy as np
import matplotlib.pyplot as plt
import os


def load_data(file_path, dtype=float):
    with open(file_path, 'r') as file:
        return np.array([[dtype(value) for value in line.split()] for line in file])


# Funkcja do tworzenia wykresu zmian S(k)(it)
def plot_functional_integral(k_values, output_dir="."):
    plt.figure(figsize=(10, 6))
    for k in k_values:
        file_path = os.path.join(output_dir, f"it_s_{k}.txt")
        data = load_data(file_path)
        iterations, s_values = data[:, 0], data[:, 1]
        plt.plot(iterations, s_values, label=f"k = {k}")
    plt.xlabel("Numer iteracji")
    plt.ylabel("Wartość całki funkcjonalnej S(k)")
    plt.title("Zmiana całki funkcjonalnej S(k) w funkcji numeru iteracji")
    plt.legend()
    plt.grid()
    plt.savefig(os.path.join(output_dir, "functional_integral_plot.png"))
    plt.show()


# Funkcja do generowania map potencjału
def plot_potential_maps(k_values, output_dir="."):
    potential_data = {}

    for k in k_values:
        file_v = os.path.join(output_dir, f"pot_{k}.txt")
        potential_data[k] = load_data(file_v)

    fig, axes = plt.subplots(1, len(k_values), figsize=(15, 5))
    for idx, k in enumerate(k_values):
        ax = axes[idx]
        im = ax.imshow(
            np.transpose(potential_data[k]),
            origin='lower',
            cmap='jet',
        )
        ax.set_title(f"k={k}")
        fig.colorbar(im, ax=ax)
        ax.set_xticks([])
        ax.set_yticks([])

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "potential_maps.png"))
    plt.show()


def main():
    output_dir = "."
    k_values = [16, 8, 4, 2, 1]

    # Tworzenie wykresów
    plot_functional_integral(k_values, output_dir)
    plot_potential_maps(k_values, output_dir)


if __name__ == "__main__":
    main()
