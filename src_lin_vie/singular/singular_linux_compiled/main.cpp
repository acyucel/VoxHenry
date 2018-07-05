#include <iostream>
#include <stdlib.h>
#include <iomanip>
#include <fstream>
#include <complex>
#include <math.h>

#include "directfn.h"

#ifndef M_PI
#define M_PI       3.14159265358979323846264338328      /* pi */
#endif

using namespace std;
///////////////////////////////////////////////////////////////////////////////

int main(int argc,char* argv[]) {

cout << "static lib" << endl;


int ker_type;
if (argc > 1)
	ker_type = atoi(argv[1]);
else ker_type = 0;
double d = 0.1;

const double r1[] = {0.0 , 0.0 , 0.0};
const double r2[] = {0.0, d, 0.0};
const double r3[] = {0.0,  d , d};
const double r4[] = {0.0 ,0.0, d};
const double r5[] = {0.0, 0.0, 2*d};
const double r6[] = {0.0, d, 2*d};


double n1[] = {-1.,0.,0.};
//double n2[] = {-1.,0.,0};
//double n3[] = {1.,0.,0.};

double rq_c[] = {d/2,d/2,0.};
double rp_c[] = {d, d/2,d/2};

	

complex <double> I_ST, I_EAo, I_EAc, I_VAo, I_VAc;
complex <double> I_ST_tri, I_EA_tri, I_VA_tri;

int N1, N2, N3, N4;

const double ko = 2*M_PI;

N1 = 10;
N2 = 5;
N3 = 5;
N4 = 5;

int size;
size = 1;



complex<double> *I = new complex<double>[size];

ofstream myfile;
    myfile.open ("Results_EAc.txt");

//create_EA(r1, r2, r3, r4, r5, r6, N1, N2, N3, N4,ko,d,rq_c,rq_c,n1,n1,ker_type, I);
//create_EA(r1,r4,r3,r2,r7,r8, N1, N2, N3, N4,ko,d,rq_c,rq_c,nq,nq,4, I_2);
//create_VA(r1,r4,r3,r4,r6,r11,r10, N1, N2, N3, N4,ko,d,rq_c,rq_c,nq,nq,4, I_3);
//I_ST = quadric_ws_st(r1, r2, r3, r4, N1, N2, N3, N4, ko,d,rq_c,rq_c,nq,np,0);
//I_EAo = quadric_ws_ea(r1,r4,r3,r2,r7,r8, N1, N2, N3, N4, ko,d,rq_c,rp_c,nq,np,ker_type);
//I_EAc = quadric_ws_ea(r1,r4,r3,r2,r5,r6, N1, N2, N3, N4, ko);

//I_VAo = quadric_ws_va(r1,r2,r3,r4,r8,r9,r10,N1, N2, N3, N4,ko);
//I_VAc = quadric_ws_va(r1,r2,r3,r4,r6,r11,r10,N1, N2, N3, N4,ko);

//I_ST_tri = directfn_ws_st(r1,r2,r3, N1, N2, N3, N4,ko);
//I_EA_tri = directfn_ws_ea(r1,r3,r2,r4, N1, N2, N3, N4,ko);
//I_VA_tri = directfn_ws_va(r3,r8,r9,r2,r1, N1, N2, N3, N4,ko);

//cout << "Runtime: "  << setprecision (4) << RunTime << " [sec]" << endl;
//cout << "I_ST = " << setprecision (20) << I_ST <<  endl ;
//cout << "I_EAo = " << setprecision (20) << I_EAo << endl ;
//for (int i = 0; i < size; i++)
//{
//	cout << "I = " << setprecision (20) << I[i]<< endl ;
//	myfile << setprecision(20) << real(I[i]) << ' ' << setprecision(20) << imag(I[i]) << endl;
//}


//cout << "I_EAc = " << setprecision (20) << I_EAc << endl ;
//cout << "I_VAo = " << setprecision (20) << I_VAo << endl ;
//cout << "I_VAc = " << setprecision (20) << I_VAc << endl ;

//cout << "I_ST_tri = " << setprecision (20) << I_ST_tri << endl ;
//cout << "I_EA_tri = " << setprecision (20) << I_EA_tri << endl ;
//cout << "I_VA_tri = " << setprecision (20) << I_VA_tri << endl ;
myfile.close();
delete I;
return 0;

}
///////////////////////////////////////////////////////////////////////////////

// End of the file

