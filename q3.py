from my_pkg.Linear import solve
from my_pkg.UTILITY import get_file
from math import log,exp

#straight line fitting 
def linear_fit(xx,yy):
    N=len(xx)  #length of data
    sum_x,sum_x2,sum_y,sum_xy,sum_y2=0,0,0,0,0 #initallising the sums
    for x,y in zip(xx,yy): #sum
        sum_x+=x
        sum_x2+=x**2
        sum_y+=y
        sum_xy+=x*y
        sum_y2+=y**2
    matrix=[[N,sum_x],[sum_x,sum_x2]] #covariance matrix
    b_vec= [sum_y,sum_xy]  #y_vector 
    params = solve(matrix,b_vec) #solving the linear equations using LU
    mean_x=sum_x/N  #mean of x
    mean_y=sum_y/N #mean of y

    #Pearson's r
    r=(sum_xy-(N*mean_x*mean_y))/(sum_x2-N*mean_x**2)**0.5/(sum_y2-N*mean_y**2)**0.5

    return params,abs(r)



mat=get_file("esem_table.dat",delimiter="\t") #output is in matrix form

x,y=[],[]
for i in mat:  
    x.append(i[0])
    y.append(i[1])

par,r=linear_fit(x,y)
w0,wc=par  #parameters are
print("Fitting with w(t)=w0+wc*t")
print("w0 = %f\nwc = %f\nr = %f\n"%(w0,wc,r))



logy=[log(j) for j in y]
par,r=linear_fit(x,logy)
w0=exp(par[0])
wc=-par[1] # minus sign is taken becase minus sign is already take in given funciton 
print("Fitting with log(w(t)) =log(w0)-wc*t")
print("w0 = %f\nwc = %f\nr = %f\n"%(w0,wc,r))

#Output
'''
Fitting with w(t)=w0+wc*t
w0 = 2.029103
wc = -0.474709
r = 0.985156

Fitting with log(w(t)) =log(w0)-wc*t
w0 = 2.204008
wc = 0.395596
r = 0.999118
'''











