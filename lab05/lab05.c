#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define DELTA 0.2
#define NX 128
#define NY 128
#define X_MAX (DELTA * NX)
#define Y_MAX (DELTA * NY)
#define TOL 1e-8

void write_potential_to_file(FILE *file, double potential[NX + 1][NY + 1], int k);

int main() {
    int k_values[] = {16, 8, 4, 2, 1};
    int k_n = sizeof(k_values) / sizeof(k_values[0]);

    double pot[NX + 1][NY + 1] = {0};

    // Ustalanie warunków brzegowych
    for (int i = 0; i <= NX; i++) {
        pot[i][NY] = -sin(2.0 * M_PI * DELTA * i / X_MAX);
        pot[i][0] = sin(2.0 * M_PI * DELTA * i / X_MAX);
    }

    for (int i = 0; i <= NY; i++) {
        pot[0][i] = sin(M_PI * DELTA * i / Y_MAX);
        pot[NX][i] = sin(M_PI * DELTA * i / Y_MAX);
    }

    int iteration = 0;

    for (int k_idx = 0; k_idx < k_n; k_idx++) {
        int k = k_values[k_idx];
        
        double prev_s = TOL;
        double next_s = 0.0;

        char fname_s[50], fname_pot[50], fname_x[50], fname_y[50];
        sprintf(fname_s, "it_s_%d.txt", k);
        sprintf(fname_pot, "pot_%d.txt", k);

        FILE *f_s = fopen(fname_s, "w");
        FILE *f_pot = fopen(fname_pot, "w");


        // Relaksacja i obliczenia dla k
        printf("Obliczanie dla k = %d\n", k);
        while (1) {
            iteration++;

            // Algorytm relaksacji
            for (int i = k; i < NX; i += k) {
                for (int j = k; j < NY; j += k) {
                    pot[i][j] = 0.25 * (pot[i + k][j] + pot[i - k][j] + pot[i][j + k] + pot[i][j - k]);
                }
            }

            next_s = 0.0;
            for (int i = 0; i <= NX-k; i += k) {
                for (int j = 0; j <= NY-k; j += k) {
                    next_s += 0.5 * k * k * DELTA * DELTA *
                        (pow((pot[i + k][j] - pot[i][j] + pot[i + k][j + k] - pot[i][j + k]) / (2.0 * k * DELTA), 2) +
                        pow((pot[i][j + k] - pot[i][j] + pot[i + k][j + k] - pot[i + k][j]) / (2.0 * k * DELTA), 2));
                }
            }

            fprintf(f_s, "%d %lf\n", iteration, next_s);
            printf("\rit %d: s = %lf", iteration, next_s);

            if (fabs((next_s - prev_s) / prev_s) < TOL) {
                break;
            }

            prev_s = next_s;
        }
        printf("\n");

        write_potential_to_file(f_pot, pot, k);

        // Zagęszczanie siatki, jeśli k != 1
        if (k != 1) {
            for (int i = 0; i < NX; i += k) {
                for (int j = 0; j < NY; j += k) {
                    pot[i + k / 2][j + k / 2] = 0.25 * (pot[i][j] + pot[i + k][j] + pot[i][j + k] + pot[i + k][j + k]);
                    pot[i + k][j + k / 2] = 0.5 * (pot[i + k][j] + pot[i + k][j + k]);
                    pot[i + k / 2][j + k] = 0.5 * (pot[i][j + k] + pot[i + k][j + k]);
                    pot[i + k / 2][j] = 0.5 * (pot[i][j] + pot[i + k][j]);
                    pot[i][j + k / 2] = 0.5 * (pot[i][j] + pot[i][j + k]);
                }
            }
        }

        for (int i = 0; i <= NX; i++) {
            pot[i][NY] = -sin(2.0 * M_PI * DELTA * i / X_MAX);
            pot[i][0] = sin(2.0 * M_PI * DELTA * i / X_MAX);
        }

        for (int i = 0; i <= NY; i++) {
            pot[0][i] = sin(M_PI * DELTA * i / Y_MAX);
            pot[NX][i] = sin(M_PI * DELTA * i / Y_MAX);
        }

        fclose(f_s);
        fclose(f_pot);
    }

    return 0;
}

// Funkcja zapisująca potencjał do pliku
void write_potential_to_file(FILE *file, double potential[NX + 1][NY + 1], int k) {
    for (int i = 0; i <= NX; i += k) {
        for (int j = 0; j <= NY; j += k) {
            fprintf(file, "%lf ", potential[i][j]);
        }
        fprintf(file, "\n");
    }
}
