from my_pkg.integrate import simpson
from math import sin,pi

L=1
g=9.8
a=sin(pi/8)  #theta_m=pi/4 so sin(theta_m/2)=sin(pi/8)

f=lambda p: 1/(1-(a*sin(p)**2))**0.5  #the function to be integrated

integration = simpson(f,0,pi/2,2,N=10)  #N is given 10 
                                        #Though fourth derviative i have given as 2
                                        #it will not be used in the simpson function
print("Time period= %f s"%(4*(L/g)**0.5*integration)) #printing the time period

#Output
#Time period= 2.256126 s
