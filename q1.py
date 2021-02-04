from my_pkg.Non_linear import newton_raphson
from math import exp
 
h=6.626e-34    #constants
k=1.381e-23
c=3e8

f=lambda x: (x-5)*exp(x)+5  #root of this function is required

def corr_b():
    error=1e-4   #initial error in x given
    x=newton_raphson(f,4.8,error)      #finding using newton raphson
    b= h*c/k/x
    for i in range(10):   #finding b upto 1e-4 precision in loop
        error/=10         
        x=newton_raphson(f,4.8,error)
        b1= h*c/k/x
        if abs(b1-b)<1e-4:  #checking wheather b has precision of 10^-4
            return b1
        b=b1
        

print("Wien's Displacement constant=%f mK"%corr_b())

#Output
#Wien's Displacement constant=0.002899 mK