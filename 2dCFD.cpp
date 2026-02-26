
#include <iostream>
#include <cmath>
#include <cstdio>
using namespace std;

double xm = 1;//領域の大きさ[m]
double zm = 1;

constexpr int zMax = 129;//格子の数
constexpr int xMax = 129;

double dx = xm / (xMax-1);//格子の大きさ[m]
double dz = zm / (zMax-1);

double keisu = 1.0/(-2.0*(1.0/(dx*dx)+1.0/(dz*dz)));//確か粘性項で使うための係数
double omega = 1.7;

double idx = 1.0 / dx;
double idz = 1.0 / dz;
double idx2 = 1.0 / (dx * dx);
double idz2 = 1.0 / (dz * dz);

static double u[zMax][xMax];
static double v[zMax][xMax];
static double p[zMax][xMax];

static double un[zMax][xMax];
static double vn[zMax][xMax];
static double pn[zMax][xMax];
//x方向の風をu
//z方向の風をv

double rho = 1;//1.225;
double nu  = 0.0001;//1.48e-5; // 空気の動粘性係数 (m^2/s)

// ファイルにヘッダ（格子情報）を書き込む
void writeHeader(FILE* fp) {
    fwrite(&zMax, sizeof(int), 1, fp);
    fwrite(&xMax, sizeof(int), 1, fp);
    fwrite(&dx, sizeof(double), 1, fp);
    fwrite(&dz, sizeof(double), 1, fp);
}

// 1フレーム分のデータ（時刻 + u,v,p 配列）を書き込む
void writeFrame(FILE* fp, double t, double u_arr[][xMax], double v_arr[][xMax], double p_arr[][xMax]) {
    fwrite(&t, sizeof(double), 1, fp);
    fwrite(u_arr, sizeof(double), zMax * xMax, fp);
    fwrite(v_arr, sizeof(double), zMax * xMax, fp);
    fwrite(p_arr, sizeof(double), zMax * xMax, fp);
}

void boundary(double u_arr[][xMax], double v_arr[][xMax], double U_lid){
    //上壁
    for(int x = 0; x < xMax; x++) {
        u_arr[0][x] = U_lid;
        v_arr[0][x] = 0;
    }
    //下壁
    for(int x = 0; x < xMax; x++) {
        u_arr[zMax-1][x] = 0;
        v_arr[zMax-1][x] = 0;
    }
    //左壁
    for(int z = 0; z < zMax; z++) {
        u_arr[z][0] = 0;
        v_arr[z][0] = 0;
    }
    //右壁
    for(int z = 0; z < zMax; z++) {
        u_arr[z][xMax-1] = 0;
        v_arr[z][xMax-1] = 0;
    }
}

int main() {
    double t_max = 50;//[s]
    double dt = 0.0005;//[s]
    double t = 0;
    double U_lid = 1;

    // 出力ファイルを開いてヘッダを書き込む
    FILE* fp = fopen("output.bin", "wb");
    if (fp == nullptr) {
        cerr << "Error: Could not open file output.bin for writing." << endl;
        return 1;
    }
    writeHeader(fp);
    double save_interval = 0.5; // 保存間隔 [s]
    double next_save = 0.0;

    int step = 0;//cout制御用

    while(t < t_max) {
        // Step 1: 仮の速度 (u*, v*) を求めるループ
        if(step % 10 == 0) {
            cout << "t:" << t << "/" << t_max << endl;
        }
        step++;
        for(int z=1;z<zMax-1;z++) {
            for(int x=1;x<xMax-1;x++) {

                //移流
                double du_dx,du_dz;
                if (u[z][x] > 0) {
                    du_dx = (u[z][x] - u[z][x-1]) * idx;
                } else {
                    du_dx = (u[z][x+1] - u[z][x]) * idx;
                }
                if (v[z][x] > 0) {
                    du_dz = (u[z][x] - u[z-1][x]) * idz;
                } else {
                    du_dz = (u[z+1][x] - u[z][x]) * idz;
                }
                double advection_u = u[z][x] * du_dx + v[z][x] * du_dz;//これがいわゆる(v・∇)v

                double dv_dx,dv_dz;
                if (u[z][x] > 0) {
                    dv_dx = (v[z][x] - v[z][x-1]) * idx;
                } else {
                    dv_dx = (v[z][x+1] - v[z][x]) * idx;
                }
                if (v[z][x] > 0) {
                    dv_dz = (v[z][x] - v[z-1][x]) * idz;
                } else {
                    dv_dz = (v[z+1][x] - v[z][x]) * idz;
                }
                double advection_v = u[z][x] * dv_dx + v[z][x] * dv_dz;
                

                //粘性項
                double du_dx2 = (u[z][x+1] - 2*u[z][x] + u[z][x-1]) * idx2;
                double du_dz2 = (u[z+1][x] - 2*u[z][x] + u[z-1][x]) * idz2;
                double viscosity_u = nu * (du_dx2 + du_dz2);

                double dv_dx2 = (v[z][x+1] - 2*v[z][x] + v[z][x-1]) * idx2;
                double dv_dz2 = (v[z+1][x] - 2*v[z][x] + v[z-1][x]) * idz2;
                double viscosity_v = nu * (dv_dx2 + dv_dz2);

                un[z][x] = u[z][x] + dt * (-advection_u + viscosity_u);
                vn[z][x] = v[z][x] + dt * (-advection_v + viscosity_v);
            }
        }
        //速度境界条件
        boundary(un,vn,U_lid);
        // Step 2: 圧力のポアソン方程式を解くループ (ヤコビ法)
        //rhs_arrを計算
        static double rhs_arr[zMax][xMax];
        //#pragma omp parallel for
        for(int z=1;z<zMax-1;z++) {
            for(int x=1;x<xMax-1;x++) {
                double du_dx = (un[z][x+1] - un[z][x-1]) * idx * 0.5;
                double dv_dz = (vn[z+1][x] - vn[z-1][x]) * idz * 0.5;
                rhs_arr[z][x] = rho/dt * (du_dx + dv_dz);
            }
        }

        int max_iterations = 3000;
        int check_interval = 10; // ★50回に1回だけ収束チェックをする

        for (int iter = 0; iter < max_iterations; ++iter) {
            double pmax = 0.0;
            bool check_convergence = (iter % check_interval == 0);

            // ==========================================
            // REDセルの更新
            // ==========================================
            for(int z=1; z<zMax-1; z++) {
                // zが偶数ならxは偶数から、zが奇数ならxは奇数からスタートし、2個飛ばしで進む
                int x_start = (z % 2 == 0) ? 2 : 1;
                for(int x=x_start; x<xMax-1; x+=2) {
                    double p_gs = keisu * (rhs_arr[z][x] - (pn[z][x+1] + pn[z][x-1]) * idx2 
                                                         - (pn[z+1][x] + pn[z-1][x]) * idz2);
                    double p_new_val = (1.0 - omega) * pn[z][x] + omega * p_gs;
                    
                    if(check_convergence) {
                        double diff = std::abs(p_new_val - pn[z][x]);
                        if(diff > pmax) pmax = diff;
                    }
                    pn[z][x] = p_new_val;
                }
            }

            // ==========================================
            // BLACKセルの更新
            // ==========================================
            for(int z=1; z<zMax-1; z++) {
                // zが偶数ならxは奇数から、zが奇数ならxは偶数からスタートし、2個飛ばしで進む
                int x_start = (z % 2 == 0) ? 1 : 2;
                for(int x=x_start; x<xMax-1; x+=2) {
                    double p_gs = keisu * (rhs_arr[z][x] - (pn[z][x+1] + pn[z][x-1]) * idx2 
                                                         - (pn[z+1][x] + pn[z-1][x]) * idz2);
                    double p_new_val = (1.0 - omega) * pn[z][x] + omega * p_gs;
                    
                    if(check_convergence) {
                        double diff = std::abs(p_new_val - pn[z][x]);
                        if(diff > pmax) pmax = diff;
                    }
                    pn[z][x] = p_new_val;
                }
            }

            if(check_convergence && pmax < 1e-4) {
                cout << "  -> : " << iter << "  (t=" << t << ")" << endl;
                break; 
            }

            // 圧力境界条件
            for(int z = 0; z < zMax; z++) {
                pn[z][0]      = pn[z][1];          // 左壁
                pn[z][xMax-1] = pn[z][xMax-2];     // 右壁
            }
            for(int x = 0; x < xMax; x++) {
                pn[0][x]      = pn[1][x];          // 上壁
                pn[zMax-1][x] = pn[zMax-2][x];     // 下壁
            }
        }
        

        // Step3:速度更新
        for(int z=1;z<zMax-1;z++) {
            for(int x=1;x<xMax-1;x++) {
                double dP_dx = (pn[z][x+1] - pn[z][x-1]) * idx * 0.5;
                double dP_dz = (pn[z+1][x] - pn[z-1][x]) * idz * 0.5;
                un[z][x] -= (dt/rho) * dP_dx;
                vn[z][x] -= (dt/rho) * dP_dz;
            }
        }

        //速度境界条件
        boundary(un,vn,U_lid);

        //gをgnで更新
        for(int z=0;z<zMax;z++) {
            for(int x=0;x<xMax;x++) {
                u[z][x] = un[z][x];
                v[z][x] = vn[z][x];
                p[z][x] = pn[z][x];
            }
        }

        t += dt;

        if(t > next_save) {
            writeFrame(fp, t, u, v,p);
            next_save += save_interval;
        }
    }
    fclose(fp);
    return 0;
}