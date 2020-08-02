#include <iostream>
#include <fstream>
#include <math.h>

using namespace std;
int main()
{
    int nper,n;
    double a[100],e[100];
    double ip,Omegap,omegap,Mp,amin,diva,emin,dive,i,Omega,omega,M;

    // *************************************** Dados de Entrada *******************************************

    // Dados do planeta em consideração, em cujo L4 se localizarão as partículas testes:
    
    // ip:	Inclinação do planeta em graus sexagesimais
    // Omegap:	Longitude do nó ascendente do planeta em graus sexagesimais
    // omegap:	Argumento do periélio do planeta em graus sexagesimais
    // Mp:	Anomalia média do planeta em graus sexagesimais

    ip=0.0; 
    Omegap=0.0;
    omegap=0.0;
    Mp=280.47;

    // Criando uma grade de partículas testes no L4 deste planeta, no espaço de parâmetros e-a, de 100x100:
    
    // amin:	Semieixo maior médio do planeta - Delta em semiexio maior
    // diva:	Divisões no semieixe maior
    // emin:	Excentricidade mínima da grade
    // dive:	Divisões na excentricidade

    amin=1.892606292092181;
    diva=0.002083368888889;
    emin=0.0;
    dive=0.5/99;
    i=ip;
    Omega=Omegap;
    omega=60+omegap;
    M=Mp;

    // ****************************************************************************************************

    // Redução de ângulo no intervalo de 0 a 2pi

    if (omega>360.0)
    {
        nper=omega/(360.0);
        omega=omega-nper*(360.0);
    }

    ofstream arquivo1("Kepler9Tesis_d_2.txt");
    n=0;
    for (int j=0;j<=99;j++)
    {
        a[j]=amin;
        for (int k=0;k<=99;k++)
        {
            n=n+1;
            e[k]=emin;
            arquivo1<<n<<"   "<<a[j]<<"   "<<e[k]<<"   "<<i<<"   "<<Omega<<"   "<<omega<<"   "<<M<< endl;
            emin=emin+dive;
        }
        amin=amin+diva;
        emin=0.0;
    }
    arquivo1.close();
    return 0;
}
