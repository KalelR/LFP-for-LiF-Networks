//====== 1 neuron
#include <iostream>
using namespace std; 

double f(double t, double y)
{
    double gl = 4.5, vr = -70.0, C = 100.0, I = 500.0;
    return (-gl*(y-vr)+I)/C;
}


double integraStep(double & y, double & t, double h, double f(double t, double y))
{
    double k1, k2, k3, k4, yn;
    k1 = h * f(t, y);
    k2 = h * f(t+h/2, y+k1/2 );
    k3 = h * f(t+h/2, y+k2/2 );
    k4 = h * f(t+h,   y+k3 );

    yn = y + (1.0/6.0)*(k1 + 2*k2 + 2*k3 + k4);
    // cout << "fun " << yn << " " << t << endl;
    y = yn;
    t += h;     
}

int main()
{
    double y, y0, t, t_step, t_exec, vr, v_th;
    int aux, aux_max, aux_trans;

    t_exec = 100;
    t_step = 0.01;   
    aux_max = t_exec/t_step;
    aux_trans = 0.0*aux_max;
    // vr = -70; v_th = 0;
    vr = 11; v_th = 18;
    y = vr+1;
    for(aux = 0; aux <= aux_max; aux++)
    {
        integraStep(y, t, t_step, f);
        if(y >= v_th){
            y = vr;
        }
        if(aux >= aux_trans)
        cout << t << " " << y << endl;
    }





    return 0;
}