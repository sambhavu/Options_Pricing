#include <iostream>
#include <vector>
#include <algorithm> 

using namespace std;

class Call {
public:
    Call(       double underlying, 
                double strike, 
                double t, 
                double iv, 
                double rate): 
                
                underlying_(underlying), 
                strike_(strike), 
                time_(t), 
                iv_(iv),
                rate_(rate){}
                
    
    virtual double get_price() {
       
       
       return 0;            //temp output
     
    }
    
private:
    double underlying_;
    double strike_;
    double time_;
    double iv_; 
    double rate_; 
};

class Put : public Call {
    
public:
    Put (
            double underlying, 
            double strike, 
            double iv, 
            double t,
            double rate
            
            ) 
    : Call  (underlying,
            strike, 
            iv, 
            t, 
            rate            ) {}

    double get_price() {
        
        return 100;
        
    }
};

int main() {
    int n;
    cin >> n;
    
    vector<Call*> options;
    
    
    for (int i = 0; i < n; ++i) {
        
        string option_type;
        
        double underlying_price; 
        double strike;
        double implied_volatility; 
        double t; 
        double i_rate; 
        
        
        cin >> option_type 
            >> underlying_price
            >> strike
            >> implied_volatility
            >> t
            >> i_rate; 
        
        if (option_type == "call") {
            options.push_back(new Call(underlying_price, strike, t, implied_volatility, i_rate));
        } 
        
        else if (option_type == "put") {
            options.push_back(new Put(underlying_price, strike, t, implied_volatility, i_rate));
        }
    }

    for (auto option : options) {
        delete option;
    }
    options.clear();

    return 0;
}
