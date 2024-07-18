#include <TMB.hpp>
// This namespace is used for printing during iterations
namespace CppAD
{
void PrintFor(const char* before, const double& var){}
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Input data
  DATA_INTEGER(amax);			// max observed age 
  DATA_INTEGER(ncrt);			// number of years
  DATA_IMATRIX(yrmat);		// matrix to reference cohort from year and age
  DATA_IVECTOR(monvec);		// matrix to reference cohort from year and age
  
  // catch data
  DATA_INTEGER(nobs);			// number of age-length observations
  DATA_IVECTOR(year);			// year of capture
  DATA_IVECTOR(MonCap);   // month of capture
  DATA_IVECTOR(cohort);		// cohort
  DATA_VECTOR(len);				// observed length
  DATA_IVECTOR(wage);			// observed whole age
  DATA_VECTOR(fage);			// fraction of month to capture (observed age = wage + fage)
  DATA_IVECTOR(gear);
  DATA_VECTOR(min_size);

    // environmental data
  DATA_INTEGER(n_beta);   // number of regression coefficients for environmental effects tested on growth increment
  DATA_INTEGER(n_alpha);  // number of regression coefficients for environmental effects tested on length-at-age-0
  DATA_MATRIX(env_dat);   // standardized environmental data
  DATA_IMATRIX(env_imat); // matrix of pointer values
  DATA_MATRIX(env_lo_dat);// environmental data for length-at-age-0
    
  // initial parameters
  // Lester growth & maturity parameters
  PARAMETER(lnmuh);		        // hyperprior mean log(h)
  PARAMETER(lnsighy);					// hyperprior sigma log(h)
  PARAMETER_VECTOR(hy);				// RE for log(h)
  PARAMETER(lnsighc);					// hyperprior sigma log(h)
  PARAMETER_VECTOR(hc);				// RE for log(h)
  PARAMETER(lnsigg);          // variation around g
  PARAMETER_VECTOR(lng);			// hyperprior mean log(g)
  PARAMETER(mulo);            // mean length-at-age-0
  PARAMETER(lnsiglo);					// variation around length-at-age-0
  PARAMETER_VECTOR(lo);				// lo - length-at-age 0
  PARAMETER(lnsig);						// log(standard deviation)
  PARAMETER(lnmua50);	        // hyperprior mean log(a50) - age @ 50% maturity
  PARAMETER(lnsiga50);				// hyperprior sigma log(a50)
  PARAMETER_VECTOR(lna50);		// RE for log(a50)
  PARAMETER_VECTOR(b_ha);     // phase 2 fixed effects
  PARAMETER_VECTOR(b_hj);     // phase 1 fixed effects
  PARAMETER_VECTOR(tcrit);    // critical temperature effect
  PARAMETER_VECTOR(alpha);    // fixed effects for l-0

  // back transformed parameters
  // Lester growth & maturity parameters
  Type sig=exp(lnsig);					    // coefficient of variation 
  vector<Type> a50=exp(lna50);	    // age-at-50% maturity

  // Derived Parameters
  matrix<Type> lenmat(ncrt,amax+1);	// storage array for length (unit, year, age)
  
  // optimization parameters
  Type nll_model=Type(0);		        // neg log likelihood of the model
  Type nll_RE=Type(0);			        // neg log likelihood of random effects
  Type nll_tot;				              // total neg log likelihood
  
  for(int crt=0;crt<ncrt;crt++)		// cohort loop
  {
    Type(lo_a)=Type(0);
    for(int j=0;j<n_alpha;j++)
    {
      lo_a+=alpha(j)*env_lo_dat(crt,j); // fixed effects for length-at-age-0
    }
    
    lenmat(crt,0)=lo(crt)+lo_a;					// length-at-age-0
    
    for(int a=1;a<=amax;a++)	// age loop
    {
      Type h=exp(hc(crt));                  // linear growth
      Type g=exp(lng(monvec(a-1)));					// g - reproductive investment 
      Type k=log(Type(1)+g/Type(3));				// von Bertalanffy k
      Type linf=Type(3)*h/exp(sum(lng)/12);	// linf - doesn't seem to work
      
      // when calculating linf using g as a month-specific covariate, we get negative growth
      // so I changed it to  calculate linf using lnmug, as this should 'integrate' reproductive
      // investment over the entire year and result in a single value for linf for each cohort. 
      
      // calculate the effects of env. covariates on growth
      // effect on adult growth increment
      // first one is same as all the rest, second is exponentiated beta to keep negative relationship
      Type(h_a)=b_ha(0)*(env_dat(env_imat(crt,a-1),0)-tcrit(0))*
        (env_dat(env_imat(crt,a-1),0)-tcrit(0));   // effect of temperature^2
      // effect on juvenile growth increment
      Type(h_j)=(b_ha(0)+b_hj(0))*(env_dat(env_imat(crt,a-1),0)-tcrit(1))*
        (env_dat(env_imat(crt,a-1),0)-tcrit(1));   // effect of temperature^2  

      for(int j=1;j<n_beta;j++)
      {  
        h_a+=b_ha(j)*env_dat(env_imat(crt,a-1),j);            // fixed effects phase 2 growth
        h_j+=(b_ha(j)+b_hj(j))*env_dat(env_imat(crt,a-1),j);  // fixed effects phase 1 growth
      }

      // calculate growth increments
      // phase 1 
      Type juv=h*exp(h_j);										// linear growth
      // individuals transitioning between age a-1 and a
      Type maturing=h*exp(h_j)*(Type(1)-(a-a50(crt)))+      // phase 1 growth component
        ((linf-(lenmat(crt,a-1)+h*exp(h_j)*(Type(1)-(a-a50(crt)))))*	// phase 2 growth component
        (Type(1)-exp(-k*(a-a50(crt)))))*
        exp(h_a);											
      // calculate phase 2 growth
      Type mature=(linf-lenmat(crt,a-1))*     // maximum growth increment 
        (Type(1)-exp(-k))*                   // proportion of max. growth increment
        exp(h_a);	              // month-specific effect on total growth increment (phase 2 growth only)
      
      // probability of being mature at age a-1
      Type pmature=Type(1)/(Type(1)+exp(-((a-Type(1))-
        a50(crt))/Type(0.1)));
      // probability of maturing between age a-1 and age a
      Type pmaturing=Type(1)/(Type(1)+
        exp(-(a-a50(crt))/Type(0.1)));
      
      // add together all growth parts * probabilities
      // TMB can't use if statements, this is a work-around
      // note that using Type(0.001, 0.01, and 0.05) did not converge properly
      // and generated odd errors
      Type inc=juv*(Type(1)-pmaturing)+	// juvenile
        maturing*(pmaturing-pmature)+		// maturing 
        mature*pmature;									// adult
      
      lenmat(crt,a)=lenmat(crt,a-1)+inc*  // maximum growth increment 
        exp(hy(yrmat(crt,a-1)));          // year-specific growth effect
    }
  }
  
  
  // catch data
  // growth in length
  for(int i=0;i<nobs;i++)
  {
    
    Type h=exp(hc(cohort(i)));              // linear growht 
    Type g=exp(lng(MonCap(i)));							// g - reproductive investment 
    Type k=log(Type(1)+g/Type(3));				  // von Bertalanffy k
    Type linf=Type(3)*h/exp(sum(lng)/12);		// linf - doesn't seem to work
    
    // when calculating linf using g as a month-specific covariate, we get negative growth
    // so I changed it to  calculate linf using lnmug, as this should 'integrate' reproductive
    // investment over the entire year and result in a single value for linf for each cohort. 
    
    
    // calculate the effects of env. covariates on growth
    // effect on adult growth increment
    Type(h_a)=b_ha(0)*(env_dat(env_imat(cohort(i),wage(i)),0)-tcrit(0))*
      (env_dat(env_imat(cohort(i),wage(i)),0)-tcrit(0));   // effect of temperature^2
    // effect on juvenile growth increment
    Type(h_j)=(b_ha(0)+b_hj(0))*(env_dat(env_imat(cohort(i),wage(i)),0)-tcrit(1))*
      (env_dat(env_imat(cohort(i),wage(i)),0)-tcrit(1));   // effect of temperature^2

    for(int j=1;j<n_beta;j++)
    {  
      h_a+=b_ha(j)*env_dat(env_imat(cohort(i),wage(i)),j);            // fixed effects phase 2
      h_j+=(b_ha(j)+b_hj(j))*env_dat(env_imat(cohort(i),wage(i)),j);  // fixed effects phase 1
    }
    
    // calculate growth increments
    // juveniles 
    Type juv=h*fage(i)*exp(h_j);										// linear growth
    // individuals maturing between age whole age and fractional age
    Type maturing=h*exp(h_j)*(Type(1)-((wage(i)+fage(i))-a50(cohort(i))))+
      ((linf-(lenmat(cohort(i),wage(i))+h*exp(h_j)*(Type(1)-((wage(i)+fage(i))-a50(cohort(i))))))*	
      (Type(1)-exp(-(k*fage(i))*((wage(i)+fage(i))-a50(cohort(i))))))*
      exp(h_a);											
    // calculate adult growth
    Type mature=(linf-lenmat(cohort(i),wage(i)))*
      (Type(1)-exp(-(k*fage(i))))*
      exp(h_a);	// von B increment 
    
    // probability of being mature at age a-1
    Type pmature=Type(1)/(Type(1)+exp(-(((wage(i)+fage(i))-Type(1))-
      a50(cohort(i)))/Type(0.1)));
    // probability of maturing bewteen age a-1 and age a
    Type pmaturing=Type(1)/(Type(1)+
      exp(-((wage(i)+fage(i))-a50(cohort(i)))/Type(0.1)));
    
    // add together all growth parts * probabilities
    // TMB can't use if statements, this is a work-around
    // note that using Type(0.001, 0.01, and 0.05) did not converge properly
    // and generated odd errors
    Type inc=juv*(Type(1)-pmaturing)+	// juvenile
      maturing*(pmaturing-pmature)+		// maturing 
      mature*pmature;	
    
    // growth increment
    Type lpred=lenmat(cohort(i),wage(i))+inc*exp(hy(year(i)));
    // Distribution not accounting for gear selectivity
    //nll_model-=dnorm(len(i),lpred,sig,TRUE);	// likelihood
    
    // Attempt at truncating distribution 
    //Type CumulativeNorm=pnorm(Type(1000),lpred,sig)-pnorm(min_size(gear(i)),lpred,sig);
    Type CumulativeNorm=Type(1)-pnorm(min_size(gear(i)),lpred,sig);
    Type TruncatedNorm=dnorm(len(i),lpred,sig,TRUE)-log(CumulativeNorm);
    
    nll_model-=TruncatedNorm;
  }
  
  // Random effects	
  // Random effects
  nll_RE-=sum(dnorm(lo,mulo,exp(lnsiglo),TRUE));	      // RE for year-effect on growth
  nll_RE-=sum(dnorm(hy,Type(0),exp(lnsighy),TRUE));	    // RE for year-effect on growth
  nll_RE-=sum(dnorm(hc,lnmuh,exp(lnsighc),TRUE));		    // RE for year-effect on growth
  nll_RE-=sum(dnorm(lna50,lnmua50,exp(lnsiga50),TRUE));	// RE for age-at-50% transition in growth curve
  // attempt at random walk for month-specific value of lng
  for(int i=1;i<12;i++)
  {
    nll_RE-=dnorm(lng(i),lng(i-1),exp(lnsigg),TRUE);
  }
  nll_RE-=dnorm(lng(0),lng(11),exp(lnsigg),TRUE);		// RE for year-effect on growth

  nll_tot=nll_model+nll_RE;		// total negative log likelihood
  
  REPORT(nll_model);	// report nll components when creating model object
  REPORT(nll_RE);	 
  REPORT(nll_tot);
  
  REPORT(lenmat);
  
  return nll_tot;		// component to minimize
}
