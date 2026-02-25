
#include <iostream>
#include <cmath>
using namespace std;

struct Cell {
    double u;
    double v;
    double p;
};

int zMax = 100;
int xMax = 100;
Cell grid[zMax][xMax];
//x方向の風をu
//z方向の風をv

double rho = 1.225;
double nu  = 1.48e-5; // 空気の動粘性係数 (m^2/s)




int main() {
    for(int z=0;z<zMax;z++) {
        for(int x=0;x<xMax;x++) {
            double u = grid[z][x].u;
            double v = grid[z][x].v;
            double p = grid[z][x].p;

            //移流
            double du_dx,du_dz;
            if u > 0 {
                du_dx = (u - grid[z][x-1].u) / dx;
            } else {
                du_dx = (grid[z][x+1].u - u) / dx;
            }

            if v > 0 {
                du_dz = (u - grid[z-1][x].u) / dz;
            } else {
                du_dz = (grid[z+1][x].u - u) / dz;
            }

            double advection_u = u * du_dx + v * du_dz;//これがいわゆる(v・∇)v

            //粘性項
            double viscosity_u = 

        }
    }
    return 0;
}