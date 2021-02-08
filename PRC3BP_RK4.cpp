/* Código que resolve o PRC3BP com RK4 */

#include <iostream>
#include <math.h>
#include <fstream>
#include <vector>

using namespace std;

/* Função que calcula as componentes da velocidade de u */
vector<double> g(vector<double> u, double t)
{
    const double mu = 0.01230312630186531;
    vector<double> v;

    v.push_back(u[2]);
    v.push_back(u[3]);
    v.push_back(2*u[3] + u[0] - (1 - mu)*(u[0] + mu)/pow(sqrt(pow(u[0], 2) + pow(u[1], 2) + pow(mu, 2) + 2*mu*u[0]), 3) - mu*(u[0] - 1 + mu)/pow(sqrt(pow(u[0], 2) + pow(u[1], 2) + pow(1 - mu, 2) - 2*(1 - mu)*u[0]), 3));
    v.push_back(-2*u[2] + u[1] - (1 - mu)*u[1]/pow(sqrt(pow(u[0], 2) + pow(u[1], 2) + pow(mu, 2) + 2*mu*u[0]), 3) - mu*u[1]/pow(sqrt(pow(u[0], 2) + pow(u[1], 2) + pow(1 - mu, 2) - 2*(1 - mu)*u[0]), 3));

    return v;
}

/* Função que calcula o produto de um escalar com um vetor */
vector<double> ProductSV(double c, vector<double> p)
{
    vector<double> cp;

    cp.push_back(p[0]*c);
    cp.push_back(p[1]*c);
    cp.push_back(p[2]*c);
    cp.push_back(p[3]*c);

    return cp;
}

/* Função que calcula a soma de dois vetores */
vector<double> Sum2V(vector<double> n, vector<double> m)
{
    vector<double> nm;

    nm.push_back(n[0] + m[0]);
    nm.push_back(n[1] + m[1]);
    nm.push_back(n[2] + m[2]);
    nm.push_back(n[3] + m[3]);

    return nm;
}

int main()
{
    int t0, i;
    double t, tf, h;
    vector<double> u;
    vector<double> k1, k2, k3, k4;
    vector<double> r1, r2, r3, r4, r5, r6, s1, s2, s3, s4, s5, s6;
    u.push_back(0.48780931);
    u.push_back(0.866025404);
    u.push_back(0);
    u.push_back(0);

    t0 = 0;
    tf = 150;
    h = 0.0001;
    i = 0;
    
    ofstream arquivo("Dados_RK4_Vetorial.txt");
    for(t = t0; t <= tf; t = t0 + i*h)
    {
    	arquivo<<t<<"   "<<u[0]<<"   "<<u[1]<<endl;

        k1 = g(u, t);

        r1 = ProductSV(h/2, k1);
        s1 = Sum2V(u, r1);
        k2 = g(s1, t + h/2);

        r2 = ProductSV(h/2, k2);
        s2 = Sum2V(u, r2);
        k3 = g(s2, t + h/2);

        r3 = ProductSV(h, k3);
        s3 = Sum2V(u, r3);
        k4 = g(s3, t + h);

        r4 = ProductSV(2.0, k2);
        r5 = ProductSV(2.0, k3);
        s4 = Sum2V(k1, r4);
        s5 = Sum2V(r5, k4);
        s6 = Sum2V(s4, s5);
        r6 = ProductSV(h/6, s6);
        u = Sum2V(u, r6);

        i = i + 1;
    }
    arquivo.close();
    
    return 0;
}
