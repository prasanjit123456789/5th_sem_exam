from random import random
from math import ceil

def call_func(f,args,*x): #helps in passing the arguments
    if args==None: return f(*x)
    else: return f(*x,args)

#-----------------*******************----------------------#

def midpoint(f,a,b,f2p_max,args=None,prec=0.001,N=None):
    if a==b: return 0
    elif a>b:
        a,b=b,a
        c=-1
    else:
        c=1
    if N==None: #if precision is given
        N=((b-a)**3 *abs(f2p_max)/24/prec)**0.5
        N=ceil(N)
        if N==0: N=1  #if linear function is given
    h=(b-a)/N
    sum=0
    for i in range(N):
        mid_x=a+i*h+h/2
        sum+=call_func(f,args,mid_x) #as weight=1 for each x_i
    return c*h*sum

#----------------------***********************-----------------------#

def trapezoidal(f,a,b,f2p_max,args=None,prec=0.001,N=None):
    if a==b: return 0
    elif a>b:
        a,b=b,a
        c=-1
    else:
        c=1
    if N==None: #if precision is given
        N=((b-a)**3 *abs(f2p_max)/12/prec)**0.5
        N=ceil(N)
        if N == 0: N = 1  # if linear function is given
    h=(b-a)/N
    sum=0
    for i in range(N+1):
        if i==0 or i==N:w=1 #weight is 1 for first and last
        else: w=2
        xi=a+i*h
        sum+=w*call_func(f,args,xi)
    return c*h*sum/2

#----------------------------*********************----------------------#

def simpson(f,a,b,f4p_max,args=None,prec=0.001,N=None):
    if a==b: return 0
    elif a>b:
        a,b=b,a
        c=-1
    else:
        c=1
    if N==None: #if precision is given
        N=((b-a)**5 *abs(f4p_max)/180/prec)**0.25
        N=ceil(N) if ceil(N)%2==0 else ceil(N)+1 # we have to divide the interval in even segments
        if N==0:N=2 #if the given function is upto 3rd degree polynomial
    h=(b-a)/N
    sum=0
    for i in range(N+1):
        if i==0 or i==N:w=1 #weight is 1 for first and last
        elif i%2==1: w=4 #for odd weight=4
        else: w=2 #for even weight=2
        xi=a+i*h
        sum+=w*call_func(f,args,xi)
    return c*h*sum/3

#---------------------------*********************--------------------#

def monte_carlo(f,ranges,args=None,N=1e6):
    #ranges contains list of ranges e.g. [(a,b),(c,d)] for a bivariate function
    sum=0
    sum2=0
    inv_pdf=1 #inverse pdf
    for tple in ranges: #for multivariate function inverse pdf is multiplication of range differences
        inv_pdf*=(tple[1]-tple[0])
    if inv_pdf==0: return 0,0
    for i in range(int(N)):
        random_set=[(tple[1]-tple[0])*random()+tple[0] for tple in ranges]
        value=call_func(f,args,*random_set)
        sum+=value
        sum2 += value ** 2
    return inv_pdf*sum/N, (sum2/N/inv_pdf - (sum/N/inv_pdf)**2)**0.5

#----------------*****************************------------------------#
