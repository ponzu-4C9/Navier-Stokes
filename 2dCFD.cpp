
#include <iostream>
#include <cmath>
using namespace std;

struct Cell {
    double u;
    double v;
    double p;
};
double xm = 10;//領域の大きさ[m]
double zm = 10;

int zMax = 100;//格子の数
int xMax = 100;

double dx = xm / (xMax-1);//格子の大きさ[m]
double dz = zm / (zMax-1);

Cell grid[zMax][xMax];
Cell grid_new[zMax][xMax];
//x方向の風をu
//z方向の風をv

double rho = 1.225;
double nu  = 1.48e-5; // 空気の動粘性係数 (m^2/s)


int main() {
    double t_max = 100;//[s]
    double dt = 0.001;//[s]
    double t = 0;

    while(t < t_max) {
        // Step 1: 仮の速度 (u*, v*) を求めるループ
        for(int z=1;z<zMax;z++) {
            for(int x=1;x<xMax;x++) {

                //移流
                double du_dx,du_dz;
                if grid[z][x].u > 0 {
                    du_dx = (grid[z][x].u - grid[z][x-1].u) / dx;
                } else {
                    du_dx = (grid[z][x+1].u - grid[z][x].u) / dx;
                }

                if grid[z][x].v > 0 {
                    dv_dz = (grid[z][x].v - grid[z-1][x].v) / dz;
                } else {
                    dv_dz = (grid[z+1][x].v - grid[z][x].v) / dz;
                }

                double advection_u = grid[z][x].u * du_dx + grid[z][x].v * du_dz;//これがいわゆる(v・∇)v
                
                double dv_dx,dv_dz;
                if grid[z][x].u > 0 {
                    dv_dx = (grid[z][x].u - grid[z][x-1].u) / dx;
                } else {
                    dv_dx = (grid[z][x+1].u - grid[z][x].u) / dx;
                }

                if grid[z][x].v > 0 {
                    dv_dz = (grid[z][x].v - grid[z-1][x].v) / dz;
                } else {
                    dv_dz = (grid[z+1][x].v - grid[z][x].v) / dz;
                }

                double advection_v = grid[z][x].u * dv_dx + grid[z][x].v * dv_dz;

                //粘性項
                double du_dx2 = (grid[z][x+1].u - 2*grid[z][x].u + grid[z][x-1].u) / (dx*dx);
                double du_dz2 = (grid[z+1][x].u - 2*grid[z][x].u + grid[z-1][x].u) / (dz*dz);
                double viscosity_u = nu * (du_dx2 + du_dz2);

                double dv_dx2 = (grid[z][x+1].v - 2*grid[z][x].v + grid[z][x-1].v) / (dx*dx);
                double dv_dz2 = (grid[z+1][x].v - 2*grid[z][x].v + grid[z-1][x].v) / (dz*dz);
                double viscosity_v = nu * (dv_dx2 + dv_dz2);

                grid_new[z][x].u = grid[z][x].u + dt * (advection_u + viscosity_u);
                grid_new[z][x].v = grid[z][x].v + dt * (advection_v + viscosity_v);
            }
        }
        // Step 2: 圧力のポアソン方程式を解くループ (ヤコビ法)
        int max_iterations = 100;
        double p_new[zMax][xMax];
        for (int iter = 0; iter < max_iterations; ++iter) {
            for(int z=1;z<zMax;z++) {
                for(int x=1;x<xMax;x++) {
                    double du_dx = (grid_new[z][x+1].u - grid_new[z][x-1].u) / (2*dx);
                    double dv_dz = (grid_new[z+1][x].v - grid_new[z-1][x].v) / (2*dz);
                    
                    double rhs = rho/dt * (du_dx + dv_dz);

                    p_new[z][x] = (grid[z][x].p + dt * rhs) / 4;
                    p_new[z][x] += (grid[z][x+1].p + grid[z][x-1].p + grid[z+1][x].p + grid[z-1][x].p) / 4;
                    p_new[z][x] *= 0.25;
                }
            }
        }
        t += dt;
    }
    return 0;
}