#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Constants
#define EPSILON 1.0
#define DELTA 0.1
#define GRID_X 150
#define GRID_Y 100
#define X_LIMIT (GRID_X * DELTA)
#define Y_LIMIT (GRID_Y * DELTA)
#define CONVERGENCE_TOLERANCE 1e-8

// Function prototypes
void initializeSource(double source[GRID_X][GRID_Y]);
double computeSourceValue(double x, double y);
void performGlobalRelaxation(const char *filename, double relaxation_factor, double source[GRID_X][GRID_Y], const char *error_file, const char *solution_file);
void performLocalRelaxation(const char *filename, double relaxation_factor, double source[GRID_X][GRID_Y]);
double calculateEnergy(double potential[GRID_X + 1][GRID_Y + 1], double source[GRID_X][GRID_Y]);
void writeGridToFile(FILE *file, double grid[GRID_X + 1][GRID_Y + 1]);

int main() {
    double source[GRID_X][GRID_Y];
    initializeSource(source);

    performGlobalRelaxation("global_0_6.txt", 0.6, source, "global_error_0_6.txt", "global_0_6_solution.txt");
    performGlobalRelaxation("global_1_0.txt", 1.0, source, "global_error_1_0.txt", "global_1_0_solution.txt");
    performLocalRelaxation("local_1_0.txt", 1.0, source);
    performLocalRelaxation("local_1_4.txt", 1.4, source);
    performLocalRelaxation("local_1_8.txt", 1.8, source);
    performLocalRelaxation("local_1_9.txt", 1.9, source);

    return 0;
}

double computeSourceValue(double x, double y) {
    double sigma_x = 0.1 * X_LIMIT;
    double sigma_y = 0.1 * Y_LIMIT;
    double peak1 = exp(-pow((x - 0.35 * X_LIMIT) / sigma_x, 2) - pow((y - 0.5 * Y_LIMIT) / sigma_y, 2));
    double peak2 = -exp(-pow((x - 0.65 * X_LIMIT) / sigma_x, 2) - pow((y - 0.5 * Y_LIMIT) / sigma_y, 2));
    return peak1 + peak2;
}

void initializeSource(double source[GRID_X][GRID_Y]) {
    for (int i = 0; i < GRID_X; i++) {
        for (int j = 0; j < GRID_Y; j++) {
            source[i][j] = computeSourceValue(i * DELTA, j * DELTA);
        }
    }
}

void performGlobalRelaxation(const char *filename, double relaxation_factor, double source[GRID_X][GRID_Y], const char *error_file, const char *solution_file) {
    FILE *file = fopen(filename, "w");
    FILE *error_file_ptr = fopen(error_file, "w");
    FILE *solution_file_ptr = fopen(solution_file, "w");

    printf("global relaxation w=%lf\n", relaxation_factor);

    double potential[GRID_X + 1][GRID_Y + 1] = {0};
    double old_potential[GRID_X + 1][GRID_Y + 1] = {0};

    for (int i = 0; i <= GRID_X; i++) {
        potential[i][0] = 10.0;
        old_potential[i][0] = 10.0;
    }

    double previous_energy = 0.0, current_energy = 0.0;
    int iteration = 0;

    while (1) {
        iteration++;

        for (int i = 1; i < GRID_X; i++) {
            for (int j = 1; j < GRID_Y; j++) {
                potential[i][j] = 0.25 * (old_potential[i + 1][j] + old_potential[i - 1][j] + old_potential[i][j + 1] +
                    old_potential[i][j - 1] + pow(DELTA, 2) / EPSILON * source[i][j]);
            }
        }

        for (int j = 1; j < GRID_Y; j++) {
            potential[0][j] = potential[1][j];
            potential[GRID_X][j] = potential[GRID_X - 1][j];
        }

        for (int i = 0; i <= GRID_X; i++) {
            for (int j = 0; j <= GRID_Y; j++) {
                old_potential[i][j] = (1.0 - relaxation_factor) * old_potential[i][j] +
                    relaxation_factor * potential[i][j];
            }
        }

        previous_energy = current_energy;
        current_energy = calculateEnergy(old_potential, source);

        fprintf(file, "%lf\n", current_energy);

        double solution_error = (previous_energy != 0.0) ? fabs((current_energy - previous_energy) / previous_energy) : 1.0;


        if (solution_error < CONVERGENCE_TOLERANCE) {
            break;
        }
    }

    writeGridToFile(solution_file_ptr, old_potential);

    double error;
    for (int i = 1; i < GRID_X; i++) {
        for (int j = 1; j < GRID_Y; j++) {
            error = (potential[i + 1][j] - 2.0 * potential[i][j] + potential[i - 1][j]) / (DELTA * DELTA) +
                (potential[i][j + 1] - 2.0 * potential[i][j] + potential[i][j - 1]) / (DELTA * DELTA) +
                source[i][j] / EPSILON;
            fprintf(error_file_ptr, "%lf ", error);
        }
        fprintf(error_file_ptr, "\n");
    }

    fclose(file);
    fclose(error_file_ptr);
    fclose(solution_file_ptr);
}

void performLocalRelaxation(const char *filename, double relaxation_factor, double source[GRID_X][GRID_Y]) {
    FILE *file = fopen(filename, "w");

    printf("local relaxation w=%lf\n", relaxation_factor);

    double potential[GRID_X + 1][GRID_Y + 1] = {0};

    for (int i = 0; i <= GRID_X; i++) {
        potential[i][0] = 10.0;
        potential[i][GRID_Y - 1] = 0.0;
    }

    double previous_energy = 0.0, current_energy = 0.0;

    while (1) {
        for (int i = 1; i < GRID_X; i++) {
            for (int j = 1; j < GRID_Y; j++) {
                potential[i][j] = (1.0 - relaxation_factor) * potential[i][j] + relaxation_factor * 0.25 * (potential[i + 1][j] +
                    potential[i - 1][j] + potential[i][j + 1] + potential[i][j - 1] + pow(DELTA, 2) / EPSILON * source[i][j]);
            }
        }

        for (int j = 1; j < GRID_Y; j++) {
            potential[0][j] = potential[1][j];
            potential[GRID_X][j] = potential[GRID_X - 1][j];
        }

        previous_energy = current_energy;
        current_energy = calculateEnergy(potential, source);

        fprintf(file, "%lf\n", current_energy);

        if (previous_energy != 0.0 &&
            fabs((current_energy - previous_energy) / previous_energy) < CONVERGENCE_TOLERANCE) {
            break;
        }
    }

    fclose(file);
}

double calculateEnergy(double potential[GRID_X + 1][GRID_Y + 1], double source[GRID_X][GRID_Y]) {
    double energy = 0.0;
    for (int i = 0; i < GRID_X; i++) {
        for (int j = 0; j < GRID_Y; j++) {
            energy += pow(DELTA, 2) *
                (0.5 * pow((potential[i + 1][j] - potential[i][j]) / DELTA, 2) + 0.5 * pow((potential[i][j + 1] -
                potential[i][j]) / DELTA, 2) - source[i][j] * potential[i][j]);
        }
    }
    return energy;
}

void writeGridToFile(FILE *file, double grid[GRID_X + 1][GRID_Y + 1]) {
    for (int j = 0; j <= GRID_Y; j++) {
        for (int i = 0; i <= GRID_X; i++) {
            fprintf(file, "%lf ", grid[i][j]);
        }
        fprintf(file, "\n");
    }
}
