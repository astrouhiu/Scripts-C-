/* Geramos um arquivo que será usado pelo swift/tools/init_tp.x, para partículas em torno 
ao ponto L4 de um determinado planeta */

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
    ip = 0; 
    Omegap = 0;
    omegap = 0;
    Mp = 280.47;  

    /* ***************************************************************************/

    int nper, n;
    double diva, dive, a, e, i, Omega, omega, M;
  
    i = ip;
    Omega = Omegap;
    omega = 60 + omegap;
    M = Mp;

    // Redução de omega no intervalo de 0 a 360 graus
    if(omega > 360)
    {
        nper = omega/(360);
        omega = omega - nper*(360);
    }

    // Passos no semieixo maior e na excentricidade da grade
    diva = 2*deltaap/99;
    dive = 0.5/99;
    
    // Semieixo maior inicial da grade	
    a = apmed - deltaap;
    n = 1;

    ofstream arquivo("Kepler-56_d_tp.txt");
    for(int j = 0; j <= 99; j++)
    {
        // Excentricidade inicial da grade
	e = 0;
        for(int k = 0; k <= 99; k++)
        {
            arquivo<<left<<setw(5)<<n<<"  "<<right<<setw(8)<<a<<"  "<<setw(12)<<e<<"  "<<setw(1)<<i<<"  "<<setw(1)<<Omega<<"  "<<setw(2)<<omega<<"  "<<setw(8)<<M<<endl;
            n = n + 1;
            e = e + dive;
        }
        a = a + diva;
    }
    arquivo.close();
    return 0;
}
