/* Gera um arquivo que será usado pelo swift/tools/init_tp.x, para partículas em 
torno ao ponto L4 de um determinado planeta */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>

using namespace std;

int main()
{
    double apmed, deltaap, ip, Omegap, omegap, Mp; 

    /* ********************************* Definir *********************************/
    /*   
    apmed:	Semieixo maior médio do planeta em au
    deltaap:	Um intervalo ao semieixo maior médio do planeta em au
    ip:	Inclinação do planeta em graus sexagesimais
    Omegap:	Longitude do nó ascendente do planeta em graus sexagesimais
    omegap:	Argumento do periélio do planeta em graus sexagesimais
    Mp:	Anomalia média do planeta em graus sexagesimais
    */
    
    apmed = 1.995743951597954;
    deltaap = 0.105143785;
    ip = 0.0; 
    Omegap = 0.0;
    omegap = 0.0;
    Mp = 280.47;  

    /* ***************************************************************************/

    int nper, n;
    double diva, dive, a, e, i, Omega, omega, M;
  
    i = ip;
    Omega = Omegap;
    omega = 60.0 + omegap;
    M = Mp;

    // Redução de omega no intervalo de 0 a 360 graus
    if(omega > 360.0)
    {
        nper = omega/360.0;
        omega = omega - nper*360.0;
    }

    // Passos no semieixo maior e na excentricidade da grade
    diva = 2.0*deltaap/99.0;
    dive = 0.5/99.0;
    
    // Semieixo maior inicial da grade	
    a = apmed - deltaap;
    n = 1;

    ofstream arquivo("Kepler-56_d_tp.txt");
    for(int j = 0; j < 100; j++)
    {
        // Excentricidade inicial da grade
	e = 0.0;
        for(int k = 0; k < 100; k++)
        {
            arquivo<<left<<setw(5)<<n<<"  "<<right<<setw(7)<<fixed<<setprecision(4)<<a<<"  "<<setw(6)<<e<<"  "<<setw(7)<<fixed<<setprecision(3)<<i<<"  "<<setw(7)<<Omega<<"  "<<setw(7)<<omega<<"  "<<setw(7)<<M<<"\n";
            n = n + 1;
            e = e + dive;
        }
        a = a + diva;
    }
    arquivo.close();
    return 0;
}
