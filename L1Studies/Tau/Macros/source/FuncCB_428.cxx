 /***************************************************************************** 
  * Project: RooFit                                                           * 
  *                                                                           * 
  * This code was autogenerated by RooClassFactory                            * 
  *****************************************************************************/ 

 // Your description goes here... 

// #include "Riostream.h" 

#include "FuncCB_428.h"
//#include "RooAbsReal.h"
//#include "/opt/exp_soft/cms/slc5_ia32_gcc434/lcg/roofit/5.25.02-cms6/include/RooAbsReal.h" 
//#include "/opt/exp_soft/cms/slc5_ia32_gcc434/lcg/roofit/5.25.02-cms6/include/RooAbsCategory.h" 
//#include "RooAbsCategory.h"
//#include "RooMath.h"

ClassImp(FuncCB) 
  
  FuncCB::FuncCB(const char *name, const char *title, 
		 RooAbsReal& _m,
		 RooAbsReal& _m0,
		 RooAbsReal& _sigma,
		 RooAbsReal& _alpha,
		 RooAbsReal& _n,
		 RooAbsReal& _norm) :
    RooAbsReal(name,title), 
    m("m","m",this,_m),
    m0("m0","m0",this,_m0),
    sigma("sigma","sigma",this,_sigma),
    alpha("alpha","alpha",this,_alpha),
    n("n","n",this,_n),
    norm("norm","norm",this,_norm)
{ 
} 


FuncCB::FuncCB(const FuncCB& other, const char* name) :  
  RooAbsReal(other,name), 
  m("m",this,other.m),
  m0("m0",this,other.m0),
  sigma("sigma",this,other.sigma),
  alpha("alpha",this,other.alpha),
  n("n",this,other.n),
  norm("norm",this,other.norm)
{ 
} 


Double_t FuncCB::valeur(Double_t et)
{
  m = et;
  return evaluate();
}

Double_t FuncCB::evaluate() const
{ 
   const double sqrtPiOver2 = 1.2533141373; // sqrt(pi/2)
   const double sqrt2 = 1.4142135624;

   Double_t sig = fabs((Double_t) sigma);
   
   Double_t t = (m - m0)/sig ;
   
   if (alpha < 0)
     t = -t;

   Double_t absAlpha = fabs(alpha / sig);
   Double_t a = TMath::Power(n/absAlpha,n)*exp(-0.5*absAlpha*absAlpha);
   Double_t b = absAlpha - n/absAlpha;

   //cout << a << " " << b << endl;

   ////// Pour la crystal ball
   // if (t <= absAlpha){
   //   return norm * exp(-0.5*t*t);
   // }
   // else
   //   {
   //     return norm * a * TMath::Power(t-b,-n) ;
   //   }

   Double_t aireGauche = (1 + ApproxErf( absAlpha / sqrt2 )) * sqrtPiOver2 ;
   Double_t aireDroite = ( a * 1/TMath::Power(absAlpha - b,n-1)) / (n - 1);
   Double_t aire = aireGauche + aireDroite;

   if ( t <= absAlpha ){
     return norm * (1 + ApproxErf( t / sqrt2 )) * sqrtPiOver2 / aire ;
   }
   else{
     return norm * (aireGauche +  a * (1/TMath::Power(t-b,n-1) - 1/TMath::Power(absAlpha - b,n-1)) / (1 - n)) / aire ;
   }
  
 } 


//_____________________________________________________________________________
Double_t FuncCB::ApproxErf(Double_t arg) const 
{
  static const double erflim = 5.0;
  if( arg > erflim )
    return 1.0;
  if( arg < -erflim )
    return -1.0;
  
  return RooMath::erf(arg);
}
