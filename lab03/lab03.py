import matplotlib.pyplot as plt


def rk2(x, v, dt, alpha):
    k1x = v
    k1v = alpha * (1 - x ** 2) * v - x

    k2x = v + dt * k1v
    k2v = alpha * (1 - (x + dt * k1x) ** 2) * (v + dt * k1v) - (x + dt * k1x)

    x_next = x + dt * 0.5 * (k1x + k2x)
    v_next = v + dt * 0.5 * (k1v + k2v)
    return x_next, v_next


def trapezoidal(x, v, dt, alpha):
    tol = 1e-10
    x_next, v_next = x, v
    while True:
        f = x_next - x - (dt / 2) * (v + v_next)
        g = v_next - v - (dt / 2) * (alpha * (1 - x ** 2) * v - x + alpha * (1 - x_next ** 2) * v_next - x_next)

        a11 = 1
        a12 = -(dt / 2)
        a21 = -(dt / 2) * (-2 * alpha * x_next * v_next - 1)
        a22 = 1 - (dt / 2) * alpha * (1 - x_next ** 2)

        det = a11 * a22 - a12 * a21
        delta_x = (a22 * -f - a12 * -g) / det
        delta_v = (a11 * -g - a21 * -f) / det

        x_next += delta_x
        v_next += delta_v

        if abs(delta_x) < tol and abs(delta_v) < tol:
            return x_next, v_next


def time_step_control(alpha, tol, method):
    x, v, dt = x0, v0, dt0
    t = 0
    t_values, x_values, v_values, dt_values = [], [], [], []
    t_values.append(t)
    x_values.append(x)
    v_values.append(v)
    dt_values.append(dt)
    while t < t_max:
        # dwa kroki dt
        x2, v2 = method(x, v, dt, alpha)
        x2, v2 = method(x2, v2, dt, alpha)

        # jeden krok 2*dt
        x1, v1 = method(x, v, 2 * dt, alpha)

        Ex = abs(x2 - x1) / (2 ** p - 1)
        Ev = abs(v2 - v1) / (2 ** p - 1)

        max_e = max(Ex, Ev)

        if max_e < tol:
            t = t + 2 * dt
            x = x2
            v = v2
            t_values.append(t)
            x_values.append(x)
            v_values.append(v)
            dt_values.append(dt)

        dt = dt * (s * tol / max_e) ** (1 / (p + 1))

    return t_values, x_values, v_values, dt_values


def simulate_and_plot(method, method_name):
    results = {}
    tolerances = [1e-2, 1e-5]

    for tol in tolerances:
        results[tol] = time_step_control(alpha, tol, method)

    plt.figure(figsize=(12, 8))

    # x(t)
    plt.subplot(2, 2, 1)
    for tol in tolerances:
        t_values, x_values, *_ = results[tol]
        plt.plot(t_values, x_values, label=f"TOL={tol}")
    plt.title(f"{method_name}: x(t)")
    plt.xlabel("t")
    plt.ylabel("x(t)")
    plt.legend()

    # v(t)
    plt.subplot(2, 2, 2)
    for tol in tolerances:
        t_values, _, v_values, _ = results[tol]
        plt.plot(t_values, v_values, label=f"TOL={tol}")
    plt.title(f"{method_name}: v(t)")
    plt.xlabel("t")
    plt.ylabel("v(t)")
    plt.legend()

    # delta_t(t)
    plt.subplot(2, 2, 3)
    for tol in tolerances:
        t_values, *_, dt_values = results[tol]
        plt.plot(t_values, dt_values, label=f"TOL={tol}")
    plt.title(f"{method_name}: delta_t(t)")
    plt.xlabel("t")
    plt.ylabel("delta_t(t)")
    plt.legend()

    # v(x)
    plt.subplot(2, 2, 4)
    for tol in tolerances:
        _, x_values, v_values, _ = results[tol]
        plt.plot(x_values, v_values, label=f"TOL={tol}")
    plt.title(f"{method_name}: v(x)")
    plt.xlabel("x")
    plt.ylabel("v(x)")
    plt.legend()

    plt.tight_layout()
    file_name = f"{method_name}_results.png"
    plt.savefig(file_name)
    plt.show()


if __name__ == "__main__":
    x0 = 0.01
    v0 = 0
    dt0 = 1
    s = 0.75
    p = 2
    t_max = 40
    alpha = 5

    simulate_and_plot(rk2, "RK2")
    simulate_and_plot(trapezoidal, "Trapezoidal")
