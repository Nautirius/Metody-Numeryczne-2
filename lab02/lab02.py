import numpy as np
import matplotlib.pyplot as plt

# Parametry
beta = 0.001
N = 500
gamma = 0.1
t_max = 100
dt = 0.1
u_0 = 1
TOL = 1e-6
max_iter = 20
alpha = beta * N - gamma


# Definicja funkcji f(t, u)
def f(u):
    return alpha * u - beta * u ** 2


# 1. Metoda trapezów
# Metoda trapezów z iteracją Picarda
def trapezoidal_picard(u0, t_max, dt, TOL, max_iter):
    t_arr = np.arange(0, t_max + dt, dt)
    u_arr = np.zeros_like(t_arr)
    u_arr[0] = u0

    for n in range(len(t_arr) - 1):
        u_prev = u_arr[n]
        u_next = u_prev  # Punkt startowy
        for _ in range(max_iter):
            u_new = u_prev + dt / 2 * (f(u_prev) + f(u_next))
            if abs(u_new - u_next) < TOL:
                u_next = u_new
                break
            u_next = u_new
        u_arr[n + 1] = u_next

    return t_arr, u_arr


# Metoda trapezów z iteracją Newtona
def trapezoidal_newton(u0, t_max, dt, TOL, max_iter):
    t_arr = np.arange(0, t_max + dt, dt)
    u_arr = np.zeros_like(t_arr)
    u_arr[0] = u0

    for n in range(len(t_arr) - 1):
        u_prev = u_arr[n]
        u_next = u_prev  # Punkt startowy
        for _ in range(max_iter):
            F = u_next - u_prev - dt / 2 * (f(u_prev) + f(u_next))
            dF = 1 - dt / 2 * (alpha - 2 * beta * u_next)
            u_new = u_next - F / dF
            if abs(u_new - u_next) < TOL:
                u_next = u_new
                break
            u_next = u_new
        u_arr[n + 1] = u_next

    return t_arr, u_arr


# 2. Niejawna metoda RK2
# Wartości z tablicy Butchera
c1 = 0.5 - np.sqrt(3) / 6
c2 = 0.5 + np.sqrt(3) / 6
a11 = 0.25
a12 = 0.25 - np.sqrt(3) / 6
a21 = 0.25 + np.sqrt(3) / 6
a22 = 0.25
b1 = 0.5
b2 = 0.5


# Metoda RK2
def implicit_rk2(u0, t_max, dt, TOL, max_iter):
    t_arr = np.arange(0, t_max + dt, dt)
    u_arr = np.zeros_like(t_arr)
    u_arr[0] = u0

    for n in range(len(t_arr) - 1):
        u_prev = u_arr[n]
        # Początkowe przybliżenie
        U1 = u_prev
        U2 = u_prev

        # Iteracja Newtona dla wyznaczenia U1 i U2
        for _ in range(max_iter):
            F1 = U1 - u_prev - dt * (a11 * f(U1) + a12 * f(U2))
            F2 = U2 - u_prev - dt * (a21 * f(U1) + a22 * f(U2))

            # Wyliczenie elementów macierzy Jacobiego
            m11 = 1 - dt * a11 * (alpha - 2 * beta * U1)
            m12 = -dt * a12 * (alpha - 2 * beta * U2)
            m21 = -dt * a21 * (alpha - 2 * beta * U1)
            m22 = 1 - dt * a22 * (alpha - 2 * beta * U2)

            # Wyznaczenie delta_U1 i delta_U2 metodą wyznacznikową
            det = m11 * m22 - m12 * m21
            delta_U1 = (F2 * m12 - F1 * m22) / det
            delta_U2 = (F1 * m21 - F2 * m11) / det

            # Aktualizacja U1 i U2
            U1_new = U1 + delta_U1
            U2_new = U2 + delta_U2

            # Sprawdzenie warunku stop
            if abs(delta_U1) < TOL and abs(delta_U2) < TOL:
                U1, U2 = U1_new, U2_new
                break

            U1, U2 = U1_new, U2_new

        # Wyznaczenie u_n+1 za pomocą wzoru korektora
        u_next = u_prev + dt * (b1 * f(U1) + b2 * f(U2))
        u_arr[n + 1] = u_next

    return t_arr, u_arr


# Funkcja do rysowania wykresów u(t) oraz z(t) = N - u(t)
def plot_results(t_arr, u_arr, title, file_name):
    z_arr = N - u_arr
    plt.figure(figsize=(8, 6))
    plt.plot(t_arr, u_arr, label="u(t) - zarażeni", color='red')
    plt.plot(t_arr, z_arr, label="z(t) - zdrowi", color='blue')
    plt.title(title)
    plt.xlabel("t")
    plt.ylabel("Liczba osób")
    plt.legend()
    plt.grid()
    plt.savefig(file_name)
    plt.show()


# Uruchamianie zadań i zapisywanie wyników
if __name__ == "__main__":
    # Zadanie 1
    # Obliczenia metodą Picarda
    t_arr, u_arr_picard = trapezoidal_picard(u_0, t_max, dt, TOL, max_iter)
    z_arr_picard = N - u_arr_picard

    plot_results(t_arr, u_arr_picard, "Metoda trapezów z iteracją Picarda", "trapezoidal_picard.png")

    # Obliczenia metodą Newtona
    _, u_arr_newton = trapezoidal_newton(u_0, t_max, dt, TOL, max_iter)
    z_arr_newton = N - u_arr_newton

    plot_results(t_arr, u_arr_newton, "Metoda trapezów z iteracją Newtona", "trapezoidal_newton.png")

    # Zadanie 2
    # Obliczenia metodą RK2
    t_arr, u_arr_rk2 = implicit_rk2(u_0, t_max, dt, TOL, max_iter)

    plot_results(t_arr, u_arr_rk2, "Niejawna metoda RK2", "implicit_rk2.png")
