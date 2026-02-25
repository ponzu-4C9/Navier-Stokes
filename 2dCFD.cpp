
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

Cell g[zMax][xMax];
Cell gn[zMax][xMax];
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
        for(int z=1;z<zMax-1;z++) {
            for(int x=1;x<xMax-1;x++) {

                //移流
                double du_dx,du_dz;
                if (g[z][x].u > 0) {
                    du_dx = (g[z][x].u - g[z][x-1].u) / dx;
                } else {
                    du_dx = (g[z][x+1].u - g[z][x].u) / dx;
                }

                if (g[z][x].v > 0) {
                    dv_dz = (g[z][x].v - g[z-1][x].v) / dz;
                } else {
                    dv_dz = (g[z+1][x].v - g[z][x].v) / dz;
                }

                double advection_u = g[z][x].u * du_dx + g[z][x].v * du_dz;//これがいわゆる(v・∇)v
                
                double dv_dx,dv_dz;
                if (g[z][x].u > 0) {
                    dv_dx = (g[z][x].u - g[z][x-1].u) / dx;
                } else {    
                    dv_dx = (g[z][x+1].u - g[z][x].u) / dx;
                }

                if (g[z][x].v > 0) {
                    dv_dz = (g[z][x].v - g[z-1][x].v) / dz;
                } else {
                    dv_dz = (g[z+1][x].v - g[z][x].v) / dz;
                }

                double advection_v = g[z][x].u * dv_dx + g[z][x].v * dv_dz;

                //粘性項
                double du_dx2 = (g[z][x+1].u - 2*g[z][x].u + g[z][x-1].u) / (dx*dx);
                double du_dz2 = (g[z+1][x].u - 2*g[z][x].u + g[z-1][x].u) / (dz*dz);
                double viscosity_u = nu * (du_dx2 + du_dz2);

                double dv_dx2 = (g[z][x+1].v - 2*g[z][x].v + g[z][x-1].v) / (dx*dx);
                double dv_dz2 = (g[z+1][x].v - 2*g[z][x].v + g[z-1][x].v) / (dz*dz);
                double viscosity_v = nu * (dv_dx2 + dv_dz2);

                gn[z][x].u = g[z][x].u + dt * (advection_u + viscosity_u);
                gn[z][x].v = g[z][x].v + dt * (advection_v + viscosity_v);
            }
        }
        // Step 2: 圧力のポアソン方程式を解くループ (ヤコビ法)
        int max_iterations = 100;
        double p_new[zMax][xMax];

        //rhs_arrを計算
        double rhs_arr[zMax][xMax];
        for(int z=1;z<zMax-1;z++) {
            for(int x=1;x<xMax-1;x++) {
                double du_dx = (gn[z][x+1].u - gn[z][x-1].u) / (2*dx);
                double dv_dz = (gn[z+1][x].v - gn[z-1][x].v) / (2*dz);
                rhs_arr[z][x] = rho/dt * (du_dx + dv_dz);
            }
        }
        //ヤコビ法で圧力のポアソン方程式を解く
        for (int iter = 0; iter < max_iterations; ++iter) {
            for(int z=1;z<zMax-1;z++) {
                for(int x=1;x<xMax-1;x++) {
                    
                    double rhs = rhs_arr[z][x];

                    double keisu = 1.0/(-2.0*(1.0/(dx*dx)+1.0/(dz*dz)));

                    p_new[z][x] = keisu * (rhs - (gn[z][x+1].p + gn[z][x-1].p)/(dx*dx) - (gn[z+1][x].p + gn[z-1][x].p)/(dz*dz));
                }
            }
            for(int z=0;z<zMax;z++) {
                for(int x=0;x<xMax;x++) {
                    gn[z][x].p = p_new[z][x];
                }
            }
        }
        //圧力境界条件
        for(int z = 0; z < zMax; z++) {
            gn[z][0].p      = gn[z][1].p;          // 左壁
            gn[z][xMax-1].p = gn[z][xMax-2].p;     // 右壁
        }
        for(int x = 0; x < xMax; x++) {
            gn[0][x].p      = gn[1][x].p;          // 上壁
            gn[zMax-1][x].p = gn[zMax-2][x].p;     // 下壁
        }

        // Step3:速度更新
        for(int z=1;z<zMax-1;z++) {
            for(int x=1;x<xMax-1;x++) {
                double dP_dx = (gn[z][x+1].p - gn[z][x-1].p) / (2*dx);
                double dP_dz = (gn[z+1][x].p - gn[z-1][x].p) / (2*dz);
                gn[z][x].u -= (dt/rho) * dP_dx;
                gn[z][x].v -= (dt/rho) * dP_dz;
            }
        }

        //速度境界条件
        double U_lid = 1.0;
        //上壁
        for(int x = 0; x < xMax; x++) {
            gn[0][x].u = U_lid;
            gn[0][x].v = 0;
        }
        //下壁
        for(int x = 0; x < xMax; x++) {
            gn[zMax-1][x].u = 0;
            gn[zMax-1][x].v = 0;
        }
        //左壁
        for(int z = 0; z < zMax; z++) {
            gn[z][0].u = 0;
            gn[z][0].v = 0;
        }
        //右壁
        for(int z = 0; z < zMax; z++) {
            gn[z][xMax-1].u = 0;
            gn[z][xMax-1].v = 0;
        }

        t += dt;
    }
    return 0;
}