#include <iostream>
#include <cmath>
#include <ctime>

using namespace std;

#define ITERATION_COUNT 8000
#define TILE_SIZE 4

int  n = 280; // количество узлов сетки
double s = 1; // размер сетки
double h = s / (double)n; // шаг
double h2 = h*h;
int fullTilesCount = n / TILE_SIZE;

inline void countValues(int lb1, int rb1, int lb2, int rb2, double* matrix)
{
    for(int j = lb1; j < rb1; ++j) {
        for(int k = lb2; k < rb2; ++k) {
            double x = h*j;
            double y = h*k;
            double f = 2*exp(x + y);
            matrix[(j * (n+1)) + k] = (matrix[((j-1) * (n+1)) + k] + matrix[((j+1) * (n+1)) + k] + matrix[((j) * (n+1)) + (k-1)] + matrix[((j) * (n+1)) + (k+1)] - h2*f) / 4;
        }
    }
}


int main()
{
    //initializition
    double *u = new double[(n+1) * (n+1)];
    for(int i = 0; i <= n+1; ++i) {
        for(int j = 0; j < n+1; ++j) {
            u[(i * (n+1)) + j] = 2.5;
        }
    }
    for(int i = 0; i <= n+1; ++i) {
        u[i] = exp(h*i);
        u[i * (n+1)] = exp(h*i);
        u[(i * (n+1)) + n] = exp(h*i + h*n);
        u[(n * (n+1)) + i] = exp(h*i + h*n);
    }

    cout << "Updated tiling Start!" << endl;
    std::clock_t clock_start = std::clock();

    //main part
    for(int i = 0; i < ITERATION_COUNT; ++i) {
        //prolog
        int j1 = 0;
        int lb1 = j1*TILE_SIZE + 1;
        int rb1 = (j1 + 1)*TILE_SIZE;
        //prolog
        int j2 = 0;
        int lb2 = j2*TILE_SIZE + 1;
        int rb2 = (j2 + 1)*TILE_SIZE;
        countValues(lb1, rb1, lb2, rb2, u);
        //main part
        for(j2 = 1; j2 < fullTilesCount; ++j2) {
            lb2 = j2*TILE_SIZE;
            rb2 = (j2 + 1)*TILE_SIZE;
            countValues(lb1, rb1, lb2, rb2, u);
        }
        //epilog
        j2 = fullTilesCount;
        lb2 = j2*TILE_SIZE;
        rb2 = n;
        countValues(lb1, rb1, lb2, rb2, u);
        //main part
        for(j1 = 1; j1 < fullTilesCount; ++j1) {
            lb1 = j1*TILE_SIZE;
            rb1 = (j1 + 1)*TILE_SIZE;
            //prolog
            j2 = 0;
            lb2 = j2*TILE_SIZE + 1;
            rb2 = (j2 + 1)*TILE_SIZE;
            countValues(lb1, rb1, lb2, rb2, u);
            //main part
            for(j2 =1; j2 < fullTilesCount; ++j2) {
                lb2 = j2*TILE_SIZE;
                rb2 = (j2 + 1)*TILE_SIZE;
                countValues(lb1, rb1, lb2, rb2, u);
            }
            //epilog
            j2 = fullTilesCount;
            lb2 = j2*TILE_SIZE;
            rb2 = n;
            countValues(lb1, rb1, lb2, rb2, u);
        }
        //epilog
        j1 = fullTilesCount;
        lb1 = j1*TILE_SIZE;
        rb1 = n;
        //prolog
        j2 = 0;
        lb2 = j2*TILE_SIZE + 1;
        rb2 = (j2 + 1)*TILE_SIZE;
        countValues(lb1, rb1, lb2, rb2, u);

        //main part
        for(j2 = 1; j2 < fullTilesCount; ++j2) {
            lb2 = j2*TILE_SIZE;
            rb2 = (j2 + 1)*TILE_SIZE;
            countValues(lb1, rb1, lb2, rb2, u);
        }
        //epilog
        j2 = fullTilesCount;
        lb2 = j2*TILE_SIZE;
        rb2 = n;
        countValues(lb1, rb1, lb2, rb2, u);
    }

    //resuts
    std::clock_t clock_finish = std::clock();
    double time = (clock_finish - clock_start) / (double)(CLOCKS_PER_SEC);
    printf("Time: %lf s\n", time);
    printf("u(%lf; %lf) = %lf\n", (n/2)*h, (n/2)*h, u[((n/2) * (n+1)) + (n/2)]);
    printf("u(%lf; %lf) = %lf\n", (n-1)*h, (n-1)*h, u[((n-1) * (n+1)) + (n-1)]);
    cout << endl;


    cout << "Only epilog Tiling start!" << endl;
    clock_start = std::clock();

    //main part
    for(int i = 0; i < ITERATION_COUNT; ++i) {
        for(int j1 = 0; j1 < fullTilesCount; ++j1) {
            int lb1 = j1*TILE_SIZE + (j1 == 0);
            int rb1 = (j1 + 1)*TILE_SIZE;
            for(int j2 =0; j2 < fullTilesCount; ++j2) {
                int lb2 = j2*TILE_SIZE + (j2 == 0);
                int rb2 = (j2 + 1)*TILE_SIZE;
                countValues(lb1, rb1, lb2, rb2, u);
            }
            //epilog
            int j2 = fullTilesCount;
            int lb2 = j2*TILE_SIZE + (j2 == 0);
            int rb2 = n;
            countValues(lb1, rb1, lb2, rb2, u);
        }
        //epilog
        int j1 = fullTilesCount;
        int lb1 = j1*TILE_SIZE + (j1 == 0);
        int rb1 = n;
        for(int j2 =0; j2 < fullTilesCount; ++j2) {
            int lb2 = j2*TILE_SIZE + (j2 == 0);
            int rb2 = (j2 + 1)*TILE_SIZE;
            countValues(lb1, rb1, lb2, rb2, u);
        }
        //epilog
        int j2 = fullTilesCount;
        int lb2 = j2*TILE_SIZE + (j2 == 0);
        int rb2 = n;
        countValues(lb1, rb1, lb2, rb2, u);
    }

    //resuts
    clock_finish = std::clock();
    time = (clock_finish - clock_start) / (double)(CLOCKS_PER_SEC);
    printf("Time: %lf s\n", time);
    printf("u(%lf; %lf) = %lf\n", (n/2)*h, (n/2)*h, u[((n/2) * (n+1)) + (n/2)]);
    printf("u(%lf; %lf) = %lf\n", (n-1)*h, (n-1)*h, u[((n-1) * (n+1)) + (n-1)]);
    cout << endl;

    cout << "Simple tiling Start!" << endl;
    clock_start = std::clock();

    //main part
    for(int i = 0; i < ITERATION_COUNT; ++i) {
        for(int j1 = 0; j1 <= fullTilesCount; ++j1) {
            int lb1 = j1*TILE_SIZE + (j1 == 0);
            int rb1 = min(n, (j1 + 1)*TILE_SIZE);
            for(int j2 =0; j2 <= fullTilesCount; ++j2) {
                int lb2 = j2*TILE_SIZE + (j2 == 0);
                int rb2 = min(n, (j2 + 1)*TILE_SIZE);
                countValues(lb1, rb1, lb2, rb2, u);
            }
        }
    }

    //resuts
    clock_finish = std::clock();
    time = (clock_finish - clock_start) / (double)(CLOCKS_PER_SEC);
    printf("Time: %lf s\n", time);
    printf("u(%lf; %lf) = %lf\n", (n/2)*h, (n/2)*h, u[((n/2) * (n+1)) + (n/2)]);
    printf("u(%lf; %lf) = %lf\n", (n-1)*h, (n-1)*h, u[((n-1) * (n+1)) + (n-1)]);
    cout << endl;


    cout << "No Tiling Start!" << endl;
    clock_start = std::clock();

    for(int i = 0; i < ITERATION_COUNT; ++i) {
        countValues(1, n, 1, n, u);
    }
    clock_finish = std::clock();
    time = (clock_finish - clock_start) / (double)(CLOCKS_PER_SEC);
    printf("Time: %lf s\n", time);

    printf("u(%lf; %lf) = %lf\n", (n/2)*h, (n/2)*h, u[((n/2) * (n+1)) + (n/2)]);
    printf("u(%lf; %lf) = %lf\n", (n-1)*h, (n-1)*h, u[((n-1) * (n+1)) + (n-1)]);

    return 0;
}

