import numpy as np
import matplotlib.pyplot as plt

# Ustawienia podstawowe
lambda_v = -1
y_0 = 1
delta_t_arr = [0.01, 0.1, 1.0]
t_lim = 5

Q_0 = 0
I_0 = 0
R = 100
L = 0.1
C = 0.001
delta_t = 0.0001
omega_0 = np.sqrt(1 / (L * C))
omega_arr = [0.5 * omega_0, 0.8 * omega_0, 1.0 * omega_0, 1.2 * omega_0]
T_0 = 2 * np.pi / omega_0
t_lim_rlc = 4 * T_0


# Funkcje do zadania 1: Rozwiązanie równania dy/dt = lambda_y

# Rozwiązanie analityczne
def analytical_method(t, lambda_v):
    return np.exp(lambda_v * t)


def get_step_num(t_lim, delta_t):
    return int(t_lim / delta_t) + 1


# Metoda Eulera
def euler_method(y_0, lambda_v, delta_t, t_lim):
    step_num = get_step_num(t_lim, delta_t)
    t_arr = np.linspace(0, t_lim, step_num)
    y_arr = np.zeros(step_num)
    y_arr[0] = y_0

    for n in range(1, step_num):
        y_arr[n] = y_arr[n - 1] + delta_t * lambda_v * y_arr[n - 1]

    return t_arr, y_arr


# Metoda RK2 (trapezów)
def rk2_method(y_0, lambda_v, delta_t, t_lim):
    step_num = get_step_num(t_lim, delta_t)
    t_arr = np.linspace(0, t_lim, step_num)
    y_arr = np.zeros(step_num)
    y_arr[0] = y_0

    for n in range(1, step_num):
        k1 = lambda_v * y_arr[n - 1]
        k2 = lambda_v * (y_arr[n - 1] + delta_t * k1)
        y_arr[n] = y_arr[n - 1] + (delta_t / 2) * (k1 + k2)

    return t_arr, y_arr


# Metoda RK4
def rk4_method(y_0, lambda_v, delta_t, t_lim):
    step_num = get_step_num(t_lim, delta_t)
    t_arr = np.linspace(0, t_lim, step_num)
    y_arr = np.zeros(step_num)
    y_arr[0] = y_0

    for n in range(1, step_num):
        k1 = lambda_v * y_arr[n - 1]
        k2 = lambda_v * (y_arr[n - 1] + (delta_t / 2) * k1)
        k3 = lambda_v * (y_arr[n - 1] + (delta_t / 2) * k2)
        k4 = lambda_v * (y_arr[n - 1] + delta_t * k3)
        y_arr[n] = y_arr[n - 1] + (delta_t / 6) * (k1 + 2 * k2 + 2 * k3 + k4)

    return t_arr, y_arr


# Generowanie wykresów dla różnych metod i kroków czasowych
def plot_methods_vs_analytical(method, method_name, t_lim):
    plt.figure(figsize=(10, 8))
    t_line = np.linspace(0, t_lim, 1000)

    for delta_t in delta_t_arr:
        t_arr, y_arr = method(y_0, lambda_v, delta_t, t_lim)
        plt.plot(t_arr, y_arr, label=f"{method_name} delta_t={delta_t}")

    y_analytical = analytical_method(t_line, lambda_v)
    plt.plot(t_line, y_analytical, label="Rozwiązanie analityczne", linestyle="--", color="black")

    plt.title(f"Porównanie {method_name} z rozwiązaniem analitycznym dla różnych kroków czasowych")
    plt.grid(True)
    plt.xlabel("Czas (t)")
    plt.ylabel("y(t)")
    plt.legend()
    plt.savefig(f"{method_name}_vs_analytical.png")
    plt.close()


# Wykres błędu globalnego dla różnych kroków czasowych
def plot_global_error(method, method_name, t_lim):
    plt.figure(figsize=(10, 8))

    for delta_t in delta_t_arr:
        t_arr, y_arr = method(y_0, lambda_v, delta_t, t_lim)
        y_analytical_line = analytical_method(t_arr, lambda_v)
        global_error = y_arr - y_analytical_line

        plt.plot(t_arr, global_error, label=f"Błąd globalny delta_t={delta_t}")

    plt.xlabel("Czas (t)")
    plt.ylabel("Błąd globalny y_num(t)-y_dok(t)")
    plt.title(f"Błąd globalny dla {method_name} dla różnych kroków czasowych")
    plt.legend()
    plt.grid(True)
    plt.savefig(f"global_error_{method_name}.png")
    plt.close()


# Funkcje do zadania 2: Rozwiązanie dla układu RLC

# Potencjał V(t) = 10 * sin(omega_v * t)
def v(t, omega_v):
    return 10 * np.sin(omega_v * t)


# Funkcja dla ładunku Q
def f_q(t, Q, I):
    return I


# Funkcja dla prądu I
def g_i(t, Q, I, R, L, C, V_t):
    return (V_t / L) - (R / L) * I - (1 / (L * C)) * Q


# Krok RK4 dla układu RLC
def rk4_step_rlc(t, Q, I, delta_t, omega_v, R, L, C):
    V_t = v(t, omega_v)
    k_Q1 = f_q(t, Q, I)
    k_I1 = g_i(t, Q, I, R, L, C, V_t)
    k_Q2 = f_q(t + delta_t / 2, Q + delta_t / 2 * k_Q1, I + delta_t / 2 * k_I1)
    k_I2 = g_i(t + delta_t / 2, Q + delta_t / 2 * k_Q1, I + delta_t / 2 * k_I1, R, L, C, v(t + delta_t / 2, omega_v))
    k_Q3 = f_q(t + delta_t / 2, Q + delta_t / 2 * k_Q2, I + delta_t / 2 * k_I2)
    k_I3 = g_i(t + delta_t / 2, Q + delta_t / 2 * k_Q2, I + delta_t / 2 * k_I2, R, L, C, v(t + delta_t / 2, omega_v))
    k_Q4 = f_q(t + delta_t, Q + delta_t * k_Q3, I + delta_t * k_I3)
    k_I4 = g_i(t + delta_t, Q + delta_t * k_Q3, I + delta_t * k_I3, R, L, C, v(t + delta_t, omega_v))
    Q_next = Q + (delta_t / 6) * (k_Q1 + 2 * k_Q2 + 2 * k_Q3 + k_Q4)
    I_next = I + (delta_t / 6) * (k_I1 + 2 * k_I2 + 2 * k_I3 + k_I4)
    return Q_next, I_next


# Rozwiązanie dla układu RLC przy różnych omega_v
def solve_rlc(omega_v, delta_t, t_lim, R, L, C, Q_0, I_0):
    t_arr = np.arange(0, t_lim, delta_t)
    Q_arr = np.zeros_like(t_arr)
    I_arr = np.zeros_like(t_arr)
    Q, I = Q_0, I_0

    for i, t in enumerate(t_arr):
        Q_arr[i] = Q
        I_arr[i] = I
        Q, I = rk4_step_rlc(t, Q, I, delta_t, omega_v, R, L, C)

    return t_arr, Q_arr, I_arr


# Generowanie wykresów dla układu RLC
def plot_rlc_solutions(omega_v_arr, delta_t, t_lim):
    plt.figure(figsize=(10, 8))

    for omega_v in omega_v_arr:
        t_arr, Q_arr, _ = solve_rlc(omega_v, delta_t, t_lim, R, L, C, Q_0, I_0)
        plt.plot(t_arr, Q_arr, label=f"omega_v={omega_v:.2f}")

    plt.title("Ładunek Q(t) dla różnych omega_v")
    plt.grid(True)
    plt.xlabel("Czas (t)")
    plt.ylabel("Ładunek Q(t)")
    plt.legend()
    plt.savefig("rlc_charge.png")
    plt.close()

    plt.figure(figsize=(10, 8))

    for omega_v in omega_v_arr:
        t_arr, _, I_arr = solve_rlc(omega_v, delta_t, t_lim, R, L, C, Q_0, I_0)
        plt.plot(t_arr, I_arr, label=f"omega_v={omega_v:.2f}")

    plt.title("Natężenie I(t) dla różnych omega_v")
    plt.grid(True)
    plt.xlabel("Czas (t)")
    plt.ylabel("Natężenie I(t)")
    plt.legend()
    plt.savefig("rlc_current.png")
    plt.close()


# Uruchamianie zadań i zapisywanie wyników
if __name__ == "__main__":
    # Wykresy dla zadania 1
    plot_methods_vs_analytical(euler_method, "Metoda Eulera", t_lim)
    plot_methods_vs_analytical(rk2_method, "Metoda RK2", t_lim)
    plot_methods_vs_analytical(rk4_method, "Metoda RK4", t_lim)

    plot_global_error(euler_method, "Metoda Eulera", t_lim)
    plot_global_error(rk2_method, "Metoda RK2", t_lim)
    plot_global_error(rk4_method, "Metoda RK4", t_lim)

    # Wykresy dla zadania 2
    plot_rlc_solutions(omega_arr, delta_t, t_lim_rlc)
