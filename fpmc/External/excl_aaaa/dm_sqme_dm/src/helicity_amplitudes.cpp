// Computes different helicity amplitudes as defined in 
// Costantini, DeTollis, Pistoni; Nuovo Cim. A2 (1971) 733-787 

#include<math.h>
#include"./functions.h"

const double PI = 4*atan(1);

void Mxxxx_fermion_dm(double x, double y, double * re, double * im){
  // some auxilliary function used in Mpppp, Mpmpm, Mpmmp.


  *re=1;
  *im=0;

  double z = - x - y;
  
  double temp;
 
  temp = 2 * ( y*y + z*z ) / ( x * x ) - 2/x;
  *re += temp * ( ReT_dm(y) + ReT_dm(z) );
  *im += temp * ( ImT_dm(y) + ImT_dm(z) );

  temp =  1/(2 * x * y) - 1/y ; 
  *re += temp * ReI_dm(x,y); 
  *im += temp * ImI_dm(x,y);

  temp =  1/(2 * x * z) - 1/z ; 
  *re += temp * ReI_dm(x,z); 
  *im += temp * ImI_dm(x,z);

  temp =  4/x +1/y +1/z + 1/(2*z*y)  -  2 * ( y*y + z*z ) / ( x * x );
  *re += temp * ReI_dm(y,z);
  *im += temp * ImI_dm(y,z);

  temp = 2* (y-z)/x;
  *re += temp * ( ReB_dm(y) - ReB_dm(z) );
  *im += temp * ( ImB_dm(y) - ImB_dm(z) );

  return;

};

void Mpppp_fermion_dm(double sred, double tred, double *re, double *im){
  // M++++ from Costantini, DeTollis, Pistoni; Nuovo Cim. A2 (1971) 733-787 

  Mxxxx_fermion_dm(sred,tred,re,im);

  return;

};

void Mpmmp_fermion_dm(double sred, double tred, double *re, double *im){
  // M+--+ from Costantini, DeTollis, Pistoni; Nuovo Cim. A2 (1971) 733-787 

  Mxxxx_fermion_dm(tred,sred,re,im);

  return;

};

void Mpmpm_fermion_dm(double sred, double tred, double *re, double *im){
// M+-+- from Costantini, DeTollis, Pistoni; Nuovo Cim. A2 (1971) 733-787 

  Mxxxx_fermion_dm(-tred-sred,tred,re,im);

  return;

};

void Mpppm_fermion_dm(double sred, double tred, double * re, double * im){
  // M+--- from Costantini, DeTollis, Pistoni; Nuovo Cim. A2 (1971) 733-787 

  double temp;

  double ured=-sred-tred;

  *re=-1;
  *im=0;

  temp=-1/sred-1/tred-1/ured;
  *re += temp*( ReT_dm(sred) + ReT_dm(tred) + ReT_dm(ured) );
  *im += temp*( ImT_dm(sred) + ImT_dm(tred) + ImT_dm(ured) );
  
  temp = 1/ured + 1/( 2 * sred * tred );
  *re += temp*ReI_dm(sred,tred);
  *im += temp*ImI_dm(sred,tred);

  temp = 1/tred + 1/( 2 * sred * ured );
  *re += temp*ReI_dm(sred,ured);
  *im += temp*ImI_dm(sred,ured);
  
  temp = 1/sred + 1/( 2 * tred * ured );
  *re += temp*ReI_dm(tred,ured);
  *im += temp*ImI_dm(tred,ured);

  return ;

};

void Mppmm_fermion_dm(double sred, double tred, double * re, double * im){
  // M++-- from Costantini, DeTollis, Pistoni; Nuovo Cim. A2 (1971) 733-787 

  double temp;
  double ured = -sred-tred;


  *re=-1;
  *im=0;
    
  temp = 1/( 2 * sred * tred );
  *re += temp*ReI_dm(sred,tred);
  *im += temp*ImI_dm(sred,tred);

  temp = 1/( 2 * sred * ured );
  *re += temp*ReI_dm(sred,ured);
  *im += temp*ImI_dm(sred,ured);
  
  temp = 1/( 2 * tred * ured );
  *re += temp*ReI_dm(tred,ured);
  *im += temp*ImI_dm(tred,ured);

  return;

};



void Mxxxx_vector_dm(double x, double y, double * re, double * im){

  // some auxilliary function used in Mpppp, Mpmpm, Mpmmp.


  *re=-1.5;
  *im=0;

  double z = - x - y;
  
  double temp;


  temp = -3* (y-z)/x;
  *re += temp * ( ReB_dm(y) - ReB_dm(z) );
  *im += temp * ( ImB_dm(y) - ImB_dm(z) );

 
  temp = -1/x*(8*x-3-6*y*z/x);
  *re += temp * ( ReT_dm(y) + ReT_dm(z) );
  *im += temp * ( ImT_dm(y) + ImT_dm(z) );

  temp =  1/x*(8*x-6-6*y*z/x)-4*(x-0.25)*(x-0.75)/(y*z);
  *re += temp * ReI_dm(y,z); 
  *im += temp * ImI_dm(y,z);


  temp = -4*(x-0.25)*(x-0.75)/(x*y); ; 
  *re += temp * ReI_dm(x,y); 
  *im += temp * ImI_dm(x,y);

  temp = -4*(x-0.25)*(x-0.75)/(x*z); ; 
  *re += temp * ReI_dm(x,z); 
  *im += temp * ImI_dm(x,z);


  return;


  
};

void Mpppp_vector_dm(double sred, double tred, double *re, double *im){

  Mxxxx_vector_dm(sred,tred,re,im);

  return;

};

void Mpmmp_vector_dm(double sred, double tred, double *re, double *im){

  Mxxxx_vector_dm(tred,sred,re,im);

  return;

};

void Mpmpm_vector_dm(double sred, double tred, double *re, double *im){

  Mxxxx_vector_dm(-tred-sred,tred,re,im);

  return;

};

void Mpppm_vector_dm(double sred, double tred, double * re, double * im){
  Mpppm_fermion_dm(sred,tred,re,im);
  *re *= -1.5;
  *im *= -1.5;

  return;
};

void Mppmm_vector_dm(double sred, double tred, double * re, double * im){

  Mppmm_fermion_dm(sred,tred,re,im);
  *re *= -1.5;
  *im *= -1.5;

  return;

};

/// Polarizable particles

//DM1 (po dim 8)

void Mpppp_DM1(double s, double t, double m, double L, double c, double *re, double *im){
double u=-s-t;
*re= -c*c/(32*PI*PI*pow(L,8))*s*s* ( ReX( s, m , L)+ ReA( t, m , L)+ ReA( u, m , L));
*im= -c*c/(32*PI*PI*pow(L,8))*s*s* ( ImX( s, m , L)+ ImA( t, m , L)+ ImA( u, m , L));
return;
};

void Mpmmp_DM1(double s, double t, double m, double L, double c, double *re, double *im){
double u=-s-t;
*re= -c*c/(32*PI*PI*pow(L,8))*t*t* ( ReX( t, m , L)+ ReA( s, m , L)+ ReA( u, m , L));
*im= -c*c/(32*PI*PI*pow(L,8))*t*t* ( ImX( t, m , L)+ ImA( s, m , L)+ ImA( u, m , L));
return;
};

void Mpmpm_DM1(double s, double t, double m, double L, double c, double *re, double *im){
double u=-s-t;
*re= -c*c/(32*PI*PI*pow(L,8))*u*u* ( ReX( u, m , L)+ ReA( t, m , L)+ ReA( s, m , L));
*im= -c*c/(32*PI*PI*pow(L,8))*u*u* ( ImX( u, m , L)+ ImA( t, m , L)+ ImA( s, m , L));
return;
};


void Mppmm_DM1(double s, double t, double m, double L, double c, double *re, double *im){
double u=-s-t;
*re= -c*c/(32*PI*PI*pow(L,8))* ( s*s*ReX( s, m , L)+ t*t*ReX( t, m , L)+ u*u*ReX( u, m , L));
*im= -c*c/(32*PI*PI*pow(L,8))* ( s*s*ImX( s, m , L)+ t*t*ImX( t, m , L)+ u*u*ImX( u, m , L));
return;
};

void Mpppm_DM1(double s, double t, double m, double L, double c, double *re, double *im){
*re= 0;
*im= 0;
return;
};

//DM2 

void Mpppp_DM2(double s, double t, double m, double L, double c, double *re, double *im){
double u=-s-t;
*re= -c*c/(8*PI*PI*pow(L,8))*s*s* ( ReC( s, m , L));
*im= -c*c/(8*PI*PI*pow(L,8))*s*s* ( ImC( s, m , L));
return;
};

void Mpmmp_DM2(double s, double t, double m, double L, double c, double *re, double *im){
double u=-s-t;
*re= -c*c/(8*PI*PI*pow(L,8))*t*t* ( ReC( t, m , L));
*im= -c*c/(8*PI*PI*pow(L,8))*t*t* ( ImC( t, m , L));
return;
};

void Mpmpm_DM2(double s, double t, double m, double L, double c, double *re, double *im){
double u=-s-t;
*re= -c*c/(8*PI*PI*pow(L,8))*u*u* ( ReC( u, m , L));
*im= -c*c/(8*PI*PI*pow(L,8))*u*u* ( ImC( u, m , L));
return;
};



void Mppmm_DM2(double s, double t, double m, double L, double c, double *re, double *im){
double u=-s-t;
*re= -c*c/(8*PI*PI*pow(L,8))* (s*s* ReC( s, m , L) + t*t* ReC( t, m , L) + u*u* ReC( u, m , L));
*im= -c*c/(8*PI*PI*pow(L,8))* (s*s* ImC( s, m , L) + t*t* ImC( t, m , L) + u*u* ImC( u, m , L));

return;
};

void Mpppm_DM2(double s, double t, double m, double L, double c, double *re, double *im){
*re= 0;
*im= 0;
return;
};


//DM3

void Mpppp_DM3(double s, double t, double m, double L, double c, double *re, double *im){
double u=-s-t;
*re= -c*c/(2*PI*PI*pow(L,4))*s*s* ( ReF0( s, m , L));
*im= -c*c/(2*PI*PI*pow(L,4))*s*s* ( ImF0( s, m , L));
return;
};

void Mpmmp_DM3(double s, double t, double m, double L, double c, double *re, double *im){
double u=-s-t;
*re= -c*c/(2*PI*PI*pow(L,4))*t*t* ( ReF0( t, m , L));
*im= -c*c/(2*PI*PI*pow(L,4))*t*t* ( ImF0( t, m , L));
return;
};

void Mpmpm_DM3(double s, double t, double m, double L, double c, double *re, double *im){
double u=-s-t;
*re= -c*c/(2*PI*PI*pow(L,4))*u*u* ( ReF0( u, m , L));
*im= -c*c/(2*PI*PI*pow(L,4))*u*u* ( ImF0( u, m , L));
return;
};

void Mppmm_DM3(double s, double t, double m, double L, double c, double *re, double *im){
double u=-s-t;
*re= -c*c/(2*PI*PI*pow(L,4))* ( s*s*ReF0( s, m , L)+ t*t*ReF0( t, m , L)+ u*u*ReF0( u, m , L));
*im= -c*c/(2*PI*PI*pow(L,4))* ( s*s*ImF0( s, m , L)+ t*t*ImF0( t, m , L)+ u*u*ImF0( u, m , L));
return;
};

void Mpppm_DM3(double s, double t, double m, double L, double c, double *re, double *im){
double u=-s-t;
*re= 0;
*im= 0;
return;
};
