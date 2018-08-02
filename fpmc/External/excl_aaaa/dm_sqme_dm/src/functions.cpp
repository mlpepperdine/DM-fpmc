#include<gsl/gsl_sf.h>
#include<math.h>

const double PI = 4*atan(1);



double ReB_dm(double z){

  double b;

  if( z==0)  // b is infinite
    {
      return 0;
    }
  else if ( z<0 )  // b > 1 
    {
      b = sqrt( 1 - 1/z );
      //     return 0.5 * b * log( (1+b) / (b-1) ) - 1;
      return b*atanh(1/b);
    }
  else if ( z>=1 ) // 0 <= b < 1 
    {
      b = sqrt( 1 - 1/z );
      //      return 0.5 * b * log( (1+b) / (1-b) ) - 1;
      return b*atanh(b);
    }
  else // b is imaginary
    {
      b = sqrt( 1/z - 1 );
      return  b * atan(1/b) ;
    }

};

double ImB_dm(double z){
  
  if(z>1)
    {
      return -PI / 2 * sqrt( 1 - 1/z );
    }
  else
    {
      return 0;
    }
};

double ReT_dm(double z){

  double b;
  double temp;

  if( z==0)  // b is infinite
    {
      return 0;
    }
  else if( z<0 ) // b>1 is real
    {
      b = sqrt( 1 - 1/z );
      temp = 0.5 * log( (1+b) / (b-1) );
      return temp*temp;
    }
  else if (z>1) // 0<b<1 is real
    {
      b = sqrt( 1 - 1/z );
      temp = 0.5 * log( (1+b) / (1-b) );
      return temp*temp - 0.25 * PI * PI;
    }
  else if (z==1) // b=0
    {
      return -0.25 * PI * PI;
    }
  else // b is imaginary
    {
      b = sqrt( 1/z - 1 );
      temp = atan( 1/b );
      return - temp*temp ;
    }

};

double ImT_dm(double z){
  
  if(z>1)
    {
      return -PI * acosh(sqrt( z ));
    }
  else
    {
      return 0;
    }
};

double Ref_dm(double q, double a){   // some auxiliary function used in ReI_dm(z,w)
                                  // (sum of the 4 dilogs)  

  double b ;
  gsl_sf_result result_re;
  gsl_sf_result result_im;

  double r;
  double theta;

  double result;

  if ( q>0 && q<1 )  // b(q) is imaginary, arguments of dilogs are complex
    {
      b = sqrt( 1/q - 1 );

      r = ( a + 1 ) / sqrt( a*a + b*b );
      theta = - atan(b/a);    // r e^itheta = (a+1)/(a+b)

      gsl_sf_complex_dilog_e(r, theta, &result_re, &result_im);
      
      result = -2 * result_re.val;  // -Re[ Li2( a+1/a+b) + Li2( a+1/a-b)] 
      
      r=(a-1) / sqrt( a * a + b * b );      // r e^itheta = (a-1)/(a+b) 
      
      gsl_sf_complex_dilog_e(r, theta, &result_re, &result_im);

      result += 2 * result_re.val;     // Re[ Li2( a-1/a+b) + Li2( a-1/a-b)] 

    }
  else // b(q) real and so are all the arguments of the dilogs
    {
      
      
      b= sqrt( 1-1/q );

      result =  gsl_sf_dilog( (a-1)/(a+b) );
      result += gsl_sf_dilog( (a-1)/(a-b) );
      result -= gsl_sf_dilog( (a+1)/(a+b) );
      result -= gsl_sf_dilog( (a+1)/(a-b) );

    }

  return result;
};

double ReI_dm(double z, double w){ // allowed regions z,w<=0 || z>=0,-z<=w<=0

  double a ;
  double lim; 

  lim = 1.0e-4;

   if ( z==0 || w==0 )
   {
      return 0;
   }
  
  // because of exact cancellations Taylor-expand 
  // the function near the origin to maintain double precision
  if (z<lim && z>-lim && w< lim && w>-lim)    
    {        
      //      return 2./3.*z*z;
      return 2./3. * z * w + 4./15. * z * w * (z+w) + 16./105. * z*w*(z*z+w*w+z*w/2);
    } 

  a = sqrt( 1-1/w-1/z);

  return ( Ref_dm(z,a) + Ref_dm(w,a) ) / (2*a);


};

double ImI_dm(double z,double w){ // allowed regions z,w<=0 || z>=0,-z<=w<=0  
  
  double b;
  double a;

  if ( z==0 || w ==0)
    {
      return 0;
    }

  a = sqrt( 1 - 1/z - 1/w );

  if ( z>1 && w<0 )
    {
      b = sqrt( 1 - 1/z );
    }
  else if ( w>1 && z<0)
    {
      b = sqrt ( 1 - 1/w );
    }
  else
    {
      return 0;  
    }

  return PI / (2*a) * log( (a-b)/(a+b) );

};


double ReF0(double  Q2, double m , double Lamb){
  double z;
  z = Q2/ (4* m*m); 
   return 2*(ReB_dm(z)-1.)+2*log(m/Lamb);
};
double ImF0(double  Q2, double m , double Lamb){
  double z;
  z = Q2/ (4* m*m); 
   return 2*ImB_dm(z);
};

double ReF1(double  Q2, double m , double Lamb){
  double z; double b2; 
    z = Q2/ (4* m*m); 
    b2 =  1 - 1/z ;
   return (0.5-b2/6.)*(ReB_dm(z)-1.)+1./18.+2./6.*log(m/Lamb);
};
double ImF1(double  Q2, double m , double Lamb){
  double z; double b2; 
    z = Q2/ (4* m*m); 
    b2 = 1 - 1/z ;
   
return (0.5-b2/6.)*ImB_dm(z);
};

double ReF2(double  Q2, double m , double Lamb){
  double z; double b2; 
    z = Q2/ (4* m*m); 
    b2 = 1 - 1/z;
   return (1./8.-b2/12.+b2*b2/40.)*(ReB_dm(z)-1.)+41./1800.-b2/120.+2./30.*log(m/Lamb);
};
double ImF2(double  Q2, double m , double Lamb){
  double z; double b2; 
    z = Q2/ (4* m*m); 
    b2 =  1 - 1/z ;
   return (1./8.-b2/12.+b2*b2/40.)*ImB_dm(z);
};
double ReA(double  Q2, double m , double Lamb){
//std::cout << "ReA" << m*m*m*m*ReF0(Q2,m,Lamb)-2*m*m*Q2*ReF1(Q2,m,Lamb)+Q2*Q2*ReF2(Q2,m,Lamb) <<"\n";
   return m*m*m*m*ReF0(Q2,m,Lamb)-2*m*m*Q2*ReF1(Q2,m,Lamb)+Q2*Q2*ReF2(Q2,m,Lamb);
};
double ImA(double  Q2, double m , double Lamb){
//std::cout << "ImA" << m*m*m*m*ImF0(Q2,m,Lamb)-2*m*m*Q2*ImF1(Q2,m,Lamb)+Q2*Q2*ImF2(Q2,m,Lamb) <<"\n";
   return m*m*m*m*ImF0(Q2,m,Lamb)-2*m*m*Q2*ImF1(Q2,m,Lamb)+Q2*Q2*ImF2(Q2,m,Lamb);
};
double ReX(double  Q2, double m , double Lamb){
//std::cout << "ReX" << (3*m*m*m*m+2*m*m*Q2)*ReF0(Q2,m,Lamb)-(30*m*m*Q2+2*Q2*Q2)*ReF1(Q2,m,Lamb)+28*Q2*Q2*ReF2(Q2,m,Lamb) <<"\n";
   return (3*m*m*m*m+2*m*m*Q2)*ReF0(Q2,m,Lamb)-(30*m*m*Q2+2*Q2*Q2)*ReF1(Q2,m,Lamb)+28*Q2*Q2*ReF2(Q2,m,Lamb);
};
double ImX(double  Q2, double m , double Lamb){
//std::cout << "ImX" << (3*m*m*m*m+2*m*m*Q2)*ImF0(Q2,m,Lamb)-(30*m*m*Q2+2*Q2*Q2)*ImF1(Q2,m,Lamb)+28*Q2*Q2*ImF2(Q2,m,Lamb) <<"\n";
   return (3*m*m*m*m+2*m*m*Q2)*ImF0(Q2,m,Lamb)-(30*m*m*Q2+2*Q2*Q2)*ImF1(Q2,m,Lamb)+28*Q2*Q2*ImF2(Q2,m,Lamb);
};
double ReC(double  Q2, double m , double Lamb){
//std::cout << "ReC" << (12*m*m*m*m+2*m*m*Q2)*ReF0(Q2,m,Lamb)-(32*m*m*Q2+2*Q2*Q2)*ReF1(Q2,m,Lamb)+24*Q2*Q2*ReF2(Q2,m,Lamb) <<"\n";
   return (12*m*m*m*m+2*m*m*Q2)*ReF0(Q2,m,Lamb)-(32*m*m*Q2+2*Q2*Q2)*ReF1(Q2,m,Lamb)+24*Q2*Q2*ReF2(Q2,m,Lamb);
};
double ImC(double  Q2, double m , double Lamb){
//std::cout << "ImC" << (12*m*m*m*m+2*m*m*Q2)*ImF0(Q2,m,Lamb)-(32*m*m*Q2+2*Q2*Q2)*ImF1(Q2,m,Lamb)+24*Q2*Q2*ImF2(Q2,m,Lamb) <<"\n";
   return (12*m*m*m*m+2*m*m*Q2)*ImF0(Q2,m,Lamb)-(32*m*m*Q2+2*Q2*Q2)*ImF1(Q2,m,Lamb)+24*Q2*Q2*ImF2(Q2,m,Lamb);
};

