from my_pkg.ode import rk4,shooter


def f(y,t):  #differential equation
    dydx = y[1]
    ddydxx = -9.8  #g=9.8m/s^2
    return [dydx,ddydxx]

#boundry values
start, end = 2, 45
#time limit
t0,t1 = 0,5
#guesses
g1,g2 = 25,35

#solving differential equation using inial velocity g1 and g2
x1,y1 = rk4(f,[start,g1],t0,t1,500)  #no of points to be taken 500 i.e dt =0.01
x2,y2 = rk4(f,[start,g2],t0,t1,500)  #no of points to be taken 500 i.e dt =0.01

h1 = y1[0][-1]  #height with first guess
h2 = y2[0][-1]  #height with second guess


#if h1 and h2 are above and below of the given height(45) then shooter method is applied
if (h1-end)*(h2-end)<0: 
    print("Heights with guesses are lying above and below the last boundry value.")
    x,y=shooter(f,start,end,g1,g2,t0,t1,500,1e-5) #shooter method g1 and g2
    print("Lunch velocity=%f m/s"%y[1][0])

#Output
'''
Heights with guesses are lying above and below the last boundry value.
Lunch velocity=33.100000 m/s
'''