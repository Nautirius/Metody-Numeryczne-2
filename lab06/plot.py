import numpy as np
import matplotlib.pyplot as plt


def load_potential_data(filename):
    data = []
    with open(filename, 'r') as file:
        for line in file:
            data.append([float(value) for value in line.split() if value])
    return np.array(data)


def plot_potential(potential, nx, ny, output_filename, title):

    plt.figure(figsize=(10, 8))
    v_min = -0.8 if 'rho' in output_filename else -10
    v_max = 0.8 if 'rho' in output_filename else 10

    plt.imshow(potential, origin='lower', cmap='seismic', extent=[0, nx, 0, ny], vmin=v_min, vmax=v_max)
    plt.colorbar(label='Potencjał (V)')
    plt.title(title)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.savefig(output_filename, dpi=300)
    plt.close()


def main():
    cases = [
        ("potential_50.txt", 50, 50, "potential_map_50.png", "Mapa potencjału (50x50)"),
        ("potential_100.txt", 100, 100, "potential_map_100.png", "Mapa potencjału (100x100)"),
        ("potential_200.txt", 200, 200, "potential_map_200.png", "Mapa potencjału (200x200)"),
        ("potential_rho_1.txt", 100, 100, "potential_map_rho_1.png", "Mapa potencjału e1=1 e2=1"),
        ("potential_rho_2.txt", 100, 100, "potential_map_rho_2.png", "Mapa potencjału e1=1 e2=2"),
        ("potential_rho_3.txt", 100, 100, "potential_map_rho_3.png", "Mapa potencjału e1=1 e2=10"),
    ]

    for filename, nx, ny, output_filename, title in cases:
        potential = load_potential_data(filename)
        plot_potential(potential, nx, ny, output_filename, title)
        print(f"Saved: {output_filename}")


if __name__ == "__main__":
    main()
