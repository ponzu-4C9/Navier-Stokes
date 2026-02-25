#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

double f(double t, double y) {
    return log(y)+exp(-t)+y*y;
}

int main() {

    double t_max = 10;
    
    scanf("%lf", &t_max);

    double dt = 0.001;
    double t = 0;


    double y = 1;
    
    double pre_t = 0;

    while(t < t_max) {
        //時間管理
        t += dt;

        //物理計算
        y += f(t,y)*dt;

        if(t - pre_t > 0.1) {
            printf("%lf→%lf\n", t, y);
            pre_t = t;
        }
    }

    return 0;
}