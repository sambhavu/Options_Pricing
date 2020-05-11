#include<algorithm>
#include<cmath>
#include<iostream>
#include<vector>
#include <conio.h>
using namespace std;

#ifndef Pi 
#define Pi 3.141592653589793238462643 
#endif 


class normal{
public: 
double normdist(double x);
double cnd(double x);
};


double normal::normdist(double x){
	double num = 0;
	double w1 = 0;
	double w2 = 0;
	
	w1 = (1/sqrt(2*Pi));
	w2 = exp(-.5*x*x);
	
	num = w1*w2; 
	
	return num;
}

double normal::cnd(double x){
	
	  double L, K, w ;
	  double const a1 = 0.31938153, 
	  a2 = -0.356563782, a3 = 1.781477937;
  double const a4 = -1.821255978, a5 = 1.330274429;

	   L = fabs(x);
	   K = 1.0 / (1.0 + 0.2316419 * L);
	   w = 1.0 - 1.0 / sqrt(2 * Pi) * exp(-L *L / 2) * (a1 * K + a2 * K *K + a3 * pow(K,3) + a4 * pow(K,4) + a5 * pow(K,5));

  if (x < 0 ){
	w= 1.0 - w;
  }
  return w;
	
}



class black_scholes{
public:
	
	normal N; 
	
	double underlying;
	double ir;
	double strike_step;
	double start_strike;
	double end_strike;
	double time;
	double v; 
	
	std::vector<double> calls;
	std::vector<double> puts; 
	
	black_scholes();
	
	void get_chain();
	
	double call(double s,double k,double sigma, double r, double t); 
	
	double put(double s,double k,double sigma, double r, double t); 
	
	double d1(double s, double k, double sigma, double r, double t); 

	void display();
};

black_scholes ::black_scholes ()
{
/*
	std::cout<<"Option Chain Generator\n\n\n"; 
	std::cout<<"Underlying Spot Price($): ";
	std::cin>> underlying;
	std::cout<<"Strike step($): ";
	std::cin>>strike_step;
	std::cout<<"Interest rate(%): ";
	std::cin>>ir;
	std::cout<<"Time(years): "; 
	std::cin>>time;
	std::cout<<"Constant Volatility(%): "; 
	std::cin>> v; 
*/
	underlying=100;
	strike_step=1;
	ir=3;
	time=.5;
	v=30;
	
	start_strike=underlying-(10*strike_step); 
	end_strike=underlying+(10*strike_step); 
	
	v=v/100;
	ir=ir/100;
	
	clrscr();
} 


void black_scholes :: display()
{
	double strike = end_strike;
 std::cout<<"PUTS\t\t\t\tstrike\t\t\t\tCALLS\n";
 int size=calls.size(); 
 
 for(int i=0; i<=size; i++) { 
 	
cout<<puts[i]<<"\t\t\t"<<strike<<"\t\t\t"<<calls[i]<<"\n"; 
 	
    if(strike==underlying)
 	{ 
 		cout<<"**********************************\n"; 
 		} 
 		
 	
 	strike-=strike_step;
 	if(strike==underlying)
 	{ 
 		cout<<"*****************ATM**************\n"; 
 		} 
 		
 }
 	
 
} 
	

double black_scholes::d1(double s, double k, double sigma, double r, double t)
{ 
	double d_1; 
	double sigma_sqr;
	double root_t;
	
	sigma_sqr=sigma*sigma;
	root_t = sqrt(t);
	
	d_1 = (log(s/k)+(r+sigma_sqr/2.0)*t)/(sigma*root_t); 
	
	return d_1; 
} 


double black_scholes :: call(double s,double k,double sigma, double r, double t)
{
	double d_1;
	double d2;
	double value;
	
	
	double root_t=sqrt(t);
	
	d_1=d1(s,k,sigma,r,t);
	
	d2=d_1-sigma*root_t; 
	
	value=s*N.normdist(d_1)-k*exp(-r*t)*N.cnd(d2);
	
	return value;
	

} 
	
	
double black_scholes :: put(double s,double k,double sigma, double r, double t)
{
	double d_1;
	double d2;
	double value;
	
	
	double root_t=sqrt(t);
	
	d_1=d1(s,k,sigma,r,t);
	
	d2=d_1-sigma*root_t; 
	
	value=k*exp(-r*t)*N.cnd(-d2)-s*N.normdist(-d_1);
	
	return value;
	

} 
	
void black_scholes :: get_chain()
{
/*
    double underlying;
	double ir;
	double strike_step;
	double start_strike;
	double end_strike;
	double time;
	double v; 
*/
	
	double call_value; 
	double put_value; 
	
	double sk=end_strike;
	
	
	while(sk>=start_strike)
	{
		call_value=call(underlying,sk,v,ir,time);
		calls.push_back(call_value);
		put_value=put(underlying,sk,v,ir,time);
		puts.push_back(put_value); 
		
		sk-=strike_step;
		
	} 
		
} 

int main()
{
	black_scholes o;
	o.get_chain();
	o.display();
	
	
	return 0;
} 
