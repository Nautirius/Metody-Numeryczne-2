#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mgmres.h"

int get_j(int l, int nx) {
    return l / (nx + 1);
}

int get_i(int l, int nx) {
    return l - get_j(l, nx) * (nx + 1);
}

void write_matrix_and_vector_to_file(const char *matrix_filename, const char *vector_filename, double *a, double *b, int nx, int N) {
    FILE *matrix_file = fopen(matrix_filename, "w");
    fprintf(matrix_file, "# l \t i_l \t j_l \t b[l]\n");
    for (int l = 0; l < N; l++) {
        fprintf(matrix_file, "%d \t %d \t %d \t %.0lf\n", l, get_i(l, nx), get_j(l, nx), b[l]);
    }
    fclose(matrix_file);

    FILE *vector_file = fopen(vector_filename, "w");
    fprintf(vector_file, "# k \t a[k]\n");
    for (int l = 0; l < N; l++) {
        fprintf(vector_file, "%d \t %.0lf\n", l, a[l]);
    }
    fclose(vector_file);
}

void write_potential_to_file(const char *potential_filename, double *V, int nx, int ny, int N) {
    FILE *output = fopen(potential_filename, "w");
    for (int i = 0; i < N - ny - 1; ++i) {
        fprintf(output, "%f ", V[i]);
        if(i % (nx + 1) == nx) { fprintf(output, "\n"); }
    }
    fclose(output);
}

double get_rho(int x, int y, double dx, double xmax, double ymax) {
    double sigma = xmax / 10.0;
    double dx_scaled = dx * x;
    double dy_scaled = dx * y;

    double rho1 = exp(-pow((dx_scaled - 0.25 * xmax) / sigma, 2) - pow((dy_scaled - 0.5 * ymax) / sigma, 2));
    double rho2 = -exp(-pow((dx_scaled - 0.75 * xmax) / sigma, 2) - pow((dy_scaled - 0.5 * ymax) / sigma, 2));
    return rho1 + rho2;
}

void solve_poisson(int nx, int ny, double dx, double e1, double e2, double v1, double v2, double v3, double v4, int rho,
                   const char *matrix_filename, const char *vector_filename, const char* potential_filename) {
    int N = (nx + 1) * (ny + 1);
    double xmax = dx * nx;
    double ymax = dx * ny;

    double *a = (double *)malloc(5 * N * sizeof(double));
    int *ia = (int *)malloc((N + 1) * sizeof(int));
    int *ja = (int *)malloc(5 * N * sizeof(int));
    double *b = (double *)malloc(N * sizeof(double));
    double *V = (double *)malloc(N * sizeof(double));

    int k = -1;
    for (int l = 0; l < N; ++l) {
        int i = get_i(l, nx);
        int j = get_j(l, nx);

        int edge = 0;
        double vb = 0.0;

        if (i == 0) { edge = 1; vb = v1; }
        if (j == ny) { edge = 1; vb = v2; }
        if (i == nx) { edge = 1; vb = v3; }
        if (j == 0) { edge = 1; vb = v4; }

        b[l] = rho ? -get_rho(get_i(l, nx), get_j(l, nx), dx, xmax, ymax) : 0.0;

        if (edge) {
            b[l] = vb;
        }

        ia[l] = -1;

        if (l - nx - 1 >= 0 && !edge) {
            k++;
            if (ia[l] < 0) ia[l] = k;
            a[k] = (i <= nx / 2) ? e1 : e2;
            a[k] /= (dx * dx);
            ja[k] = l - nx - 1;
        }

        if (l - 1 >= 0 && !edge) {
            k++;
            if (ia[l] < 0) ia[l] = k;
            a[k] = (i <= nx / 2) ? e1 : e2;
            a[k] /= (dx * dx);
            ja[k] = l - 1;
        }

        k++;
        if (ia[l] < 0) ia[l] = k;
        if (!edge) {
            double eps = (i <= nx / 2) ? e1 : e2;
            double eps_r = get_i(l + 1, nx) <= nx / 2 ? e1 : e2;
            double eps_u = get_i(l + nx + 1, nx) <= ny / 2 ? e1 : e2;
            a[k] = -(2 * eps + eps_r + eps_u) / (dx * dx);
        } else {
            a[k] = 1.0;
        }
        ja[k] = l;

        if (l < N && !edge) {
            k++;
            a[k] = (get_i(l + 1, nx) <= nx / 2) ? e1 : e2;
            a[k] /= (dx * dx);
            ja[k] = l + 1;
        }

        if (l < N - nx - 1 < N && !edge) {
            k++;
            a[k] = (get_i(l + nx + 1, nx) <= nx / 2) ? e1 : e2;
            a[k] /= (dx * dx);
            ja[k] = l + nx + 1;
        }
    }

    ia[N] = k + 1;

    pmgmres_ilu_cr(N, k + 1, ia, ja, a, V, b, 500, 500, 1e-8, 1e-8);

    write_matrix_and_vector_to_file(matrix_filename, vector_filename, a, b, nx, N);

    write_potential_to_file(potential_filename, V, nx, ny, N);

    free(a);
    free(ia);
    free(ja);
    free(b);
    free(V);
}

int main() {
    solve_poisson(4, 4, 0.1, 1, 1, 10, -10, 10, -10, 0, "matrix_4.txt", "vector_4.txt", "potential_4.txt");

    solve_poisson(50, 50, 0.1, 1, 1, 10, -10, 10, -10, 0, "matrix_50.txt", "vector_50.txt", "potential_50.txt");
    solve_poisson(100, 100, 0.1, 1, 1, 10, -10, 10, -10, 0, "matrix_100.txt", "vector_100.txt", "potential_100.txt");
    solve_poisson(200, 200, 0.1, 1, 1, 10, -10, 10, -10, 0, "matrix_200.txt", "vector_200.txt", "potential_200.txt");

    // rho
    solve_poisson(100, 100, 0.1, 1, 1, 0, 0, 0, 0, 1, "matrix_rho_1.txt", "vector_rho_1.txt", "potential_rho_1.txt");
    solve_poisson(100, 100, 0.1, 1, 2, 0, 0, 0, 0, 1, "matrix_rho_2.txt", "vector_rho_2.txt", "potential_rho_2.txt");
    solve_poisson(100, 100, 0.1, 1, 10, 0, 0, 0, 0, 1, "matrix_rho_3.txt", "vector_rho_3.txt", "potential_rho_3.txt");
    return 0;
}
