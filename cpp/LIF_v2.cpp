//======  uncoupled network
#include <iostream>
#include <vector>
using namespace std; 

double f(double t, double y)
{
    double gl = 4.5, vr = -70.0, C = 100.0, I = 500.0;
    return (-gl*(y-vr)+I)/C;
}


double integraStep(vector<double> & v_y, double & t, double h, double f(double t, double y))
{
    double k1, k2, k3, k4, y;
    int N = 2;
    for(int i = 0; i < N; i++)
    {
        y = v_y[i];
        k1 = h * f(t, y);
        k2 = h * f(t+h/2, y+k1/2 );
        k3 = h * f(t+h/2, y+k2/2 );
        k4 = h * f(t+h,   y+k3 );

        y = y + (1.0/6.0)*(k1 + 2*k2 + 2*k3 + k4);
        // cout << "fun " << yn << " " << t << endl;
        // y = yn;
        v_y[i] = y;
        // cout << v_y[i] << endl;
    }
        t += h;     
}

int main()
{
    double y, y0, t, t_step, t_exec, vr, v_th;
    int aux, aux_max, aux_trans, N, i;

    t_exec = 100;
    t_step = 0.01;   
    t = 0;
    aux_max = t_exec/t_step;
    aux_trans = 0.0*aux_max;
    // vr = -70; v_th = 0;
    vr = 11; v_th = 18;
    // y = vr+1;
    N = 1000;
    vector<double> v_y(N, vr+1);
    v_y[1] = vr+5;
    for(aux = 0; aux <= aux_max; aux++)
    {
        integraStep(v_y, t, t_step, f);
        for(i = 0; i < N; i++){
            if(v_y[i] >= v_th){
                v_y[i] = vr;
            }
        }
     
        // cout << t << " ";
        // cout << v_y[i] << " ";
        // cout << endl;
    }





    return 0;
}