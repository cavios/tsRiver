// Separable covariance on 2D lattice with AR1 structure in each direction.
#include <TMB.hpp>

/* Parameter transform */
template <class Type>
vector <Type> f(vector<Type> x){return Type(2)/(Type(1) + exp(-Type(2) * x)) - Type(1);}

template <class Type>
Type ilogit(Type x){
  return Type(1.0)/(Type(1.0)+exp(-x));
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(H)
  DATA_VECTOR(x)
  DATA_VECTOR(w)    
  DATA_VECTOR(xknot)
  DATA_VECTOR(VS)
  DATA_IVECTOR(tidx)
  DATA_IVECTOR(satid)
   

  PARAMETER_VECTOR(eta);
  PARAMETER_VECTOR(mu);
  PARAMETER_VECTOR(logSd);
  PARAMETER_VECTOR(logAscale);
  PARAMETER_VECTOR(bias);
  PARAMETER_VECTOR(transf_phi);
  vector<Type> phi=f(transf_phi);
  ADREPORT(phi);
  vector<Type> sd=exp(logSd);
  
  vector<Type> cumExpMu(mu.size());
  cumExpMu(0)=mu(0); 
  
  for(int i=1; i<mu.size(); ++i)cumExpMu(i)=cumExpMu(i-1)+exp(mu(i));
  using namespace density;
  using namespace tmbutils;
  Type res=0;  

  //process part
  res += -dnorm(eta(0),Type(0),Type(0.00001),true);
  for(int i=1;i<eta.size();i++){    
    res += -dnorm(eta(i),phi(0)*eta(i-1),sqrt(Type(1)/(1-phi(0)*phi(0))),true);
  }
  
  splinefun<Type> sp(xknot,cumExpMu);
  splinefun<Type> spS(xknot,exp(logAscale));

  vector<Type> pred(H.size());
  vector<Type> predsd(H.size());
  vector<Type> residual(H.size());

  //observation likelihood
  for(int i=0;i<H.size();i++){
    
    pred(i)=eta(tidx(i)-1)*spS(x(i))+sp(x(i));
    if(satid(i)<bias.size()){
      pred(i)+=bias(satid(i));
    }    
  }

  // weighted obs likelihhod
  for(int i=0;i<H.size();i++){
    predsd(i)=sd(satid(i));
    res+= -w(i)*dnorm(H(i),pred(i),predsd(i),true); 
  }
  //residuals
  for(int i=0;i<H.size();i++)residual(i)=(H(i)-pred(i))/sd(satid(i));
  REPORT(pred);  REPORT(predsd);
  REPORT(cumExpMu);
  REPORT(residual);
  REPORT(eta);
  

  //calculate water level ts at VS
  vector<Type> hVS(eta.size());
  for (int j=0;j<eta.size();j++){
    hVS(j)=eta(j)*spS(VS(0))+sp(VS(0));
  }
  //get uncertainty of hVS
  ADREPORT(hVS);
  return res;
}
