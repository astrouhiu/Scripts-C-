/* This code generates the "big.in" file */

#include <iostream>
#include <math.h>
#include <fstream>
#include <iomanip>

using namespace std;

/* Function that converts kg to solar masses */
double kg_to_solar_mass(double mass)
{
    return mass/(1.988500*pow(10, 30));
}

int main()
{
    /* ********************************** Define **********************************/
    /*
    num:	Number of big bodies, except the star
    name:	Names of bodies, each having up to 8 characters
    m:    	Masses of bodies in solar masses
    r:     	Maximum distances from the bodies (in Hill radii) that constitute a 
    		close encounter
    d:     	Densities of bodies in g/cm**3
    a:     	Semimajor axes in au
    e:     	Eccentricities
    I:     	Inclinations in degrees
    g:     	Arguments of pericenter in degrees
    n:     	Longitudes of the ascending node in degrees
    M:     	Mean anomalies in degrees
    */

    int num = 4;
    string name[] = {"Jupiter", "Saturn", "Uranus", "Neptune"};
    double m[] = {kg_to_solar_mass(1898.13*pow(10, 24)), kg_to_solar_mass(5.68319*pow(10, 26)), kg_to_solar_mass(86.8103*pow(10, 24)), kg_to_solar_mass(102.41*pow(10, 24))};
    double r[] = {3, 3, 3, 3};
    double d [] = {1.326, 0.687, 1.318, 1.638};
    // In this case, the orbital elements are given with respect to the ecliptic
    double a[] = {5.200755881510799, 9.572725520820169, 1.913070406983577*pow(10, 1), 3.002018230375301*pow(10, 1)};
    double e[] = {4.917281465052510*pow(10, -2), 5.150805026849767*pow(10, -2), 4.958554320329014*pow(10, -2), 6.342637279537105*pow(10, -3)};
    double I[] = {1.303686643355775, 2.489894958601465, 7.704542685081738*pow(10, -1), 1.778510161870306};
    double g[] = {2.737878544953696*pow(10, 2), 3.402348185337537*pow(10, 2), 9.903724553899913*pow(10, 1), 2.783731686579511*pow(10, 2)};
    double n[] = {1.005204718236170*pow(10, 2), 1.135533964414320*pow(10, 2), 7.410487245594420*pow(10, 1), 1.319381011067463*pow(10, 2)};
    double M[] = {1.927364306651317*pow(10, 2), 1.706166289887654*pow(10, 2), 2.152650431101928*pow(10, 2), 2.926977774060805*pow(10, 2)};

    /* ****************************************************************************/

    string epoch;
    string str1, str2, str3, str4, str5;
    
    epoch = "0";
    
    str1 = ")O+_06 Big-body initial data  (WARNING: Do not delete this line!!)";
    str2 = ") Lines beginning with `)' are ignored.";
    str3 = ")---------------------------------------------------------------------";
    str4 = " style (Cartesian, Asteroidal, Cometary) = Ast";
    str5 = " epoch (in days) = " + epoch;

    ofstream arquivo("big.in");
    arquivo<<str1<<"\n"<<str2<<"\n"<<str3<<"\n"<<str4<<"\n"<<str5<<"\n"<<str3<<"\n";
    for(int i = 0; i < num; i++)
    {
        arquivo<<setw(9)<<name[i]<<"  "<<"m="<<setprecision(8)<<scientific<<m[i]<<"  "<<"r="<<fixed<<setprecision(1)<<r[i]<<"  "<<"d="<<fixed<<setprecision(3)<<d[i]<<"\n"<<setw(14)<<fixed<<setprecision(8)<<a[i]<<"  "<<setw(10)<<e[i]<<"  "<<setw(8)<<fixed<<setprecision(4)<<I[i]<<"  "<<setw(8)<<g[i]<<"  "<<setw(8)<<n[i]<<"  "<<setw(8)<<M[i]<<"  0  0  0\n";
    }
    arquivo.close();	
    return 0;
}
