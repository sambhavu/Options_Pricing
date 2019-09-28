#include <cmath> 
#include <algorithm>
#include <iostream>
#include <string>
using namespace std; 


class option{
    public: 
    double MertonJump( 
                    double S,
                    double X, 
                    double r, 
                    double sigma, 
                    double time_to_maturity,
                    double lambda, 
                    double kappa, 
                    double delta);
    
}



double option :: MertonJump( 
                    double S,
                    double X, 
                    double r, 
                    double sigma, 
                    double time_to_maturity,
                    double lambda, 
                    double kappa, 
                    double delta)
                    
{ 
    const int maxn = 50; 
    double tau = time_to_maturity;
    double sigma_sqr = sigma*sigma      
    double delta_sqr = delta*delta;
    double lambdaprime = lambda * (1+kappa); 
    double gamma = log(1+kappa);
    double c = exp(-lambdaprime*tau)*option_price_call(S,X,r-lambda*kappa,sigma,tau);
    double log_n = 0;
    
    for(int n = 1; n<=maxn; ++n){
        log_n+=log(double(n)); 
        double sigma_n = sqrt(sigma_sqr+n*delta_sqr/tau;
        double r_n = r-lambda*kappa+n*gamma/tau;
        c += exp(-lambdaprime*tau+n*log(lambdaprime*tau)-log_n)*
        option_price_call(S,X,r_n,sigma_n,tau);
    }
    
    return c;    
}






int main()
{
    double S = 100; 
    double K = 100; 
    double r = .05; 
    double sigma = .3; 
    double time_to_maturity = 1;
    double lambda = .5;
    double kappa = .5; 
    double delta = .5; 
    
    option underlying; 
    underlying.MertonJump(S,K,r,sigma,time_to_maturity,lambda,kappa,delta);
    return 0;
    
}
