from .Non_linear import newton_raphson,call_func

def mul(a,lst): #multiplication of elements of a list with a scalar
    return [a*i for i in lst]
def add(*lst): #list addition similar to vector addition
    return [sum(i) for i in zip(*lst)]


#---------------*********************--------------------#

def euler_forward(f,y0,x0,xstop,N,args=None):
    h = (xstop - x0) / N
    result = [y0];c=[x0]
    iter=0
    while iter<N:
        iter+=1
        d = call_func(f, args, result[-1], x0)
        k1 = h*d
        yn =result[-1]+ k1
        result.append(yn)
        x0+=h
        c.append(x0)
    return c,result

#---------------------******************---------------------#

##The function to get the y(x+h) in short it converts the differential form to nonlinear form
##Then it will be fed to Newton_Raphson
def func(t,args):
    # args=[f,y0,x,h,arg]
    c=call_func(args[0],args[4],t,args[2]+args[3])
    return -t+args[1]+args[3]*c

#--------------------**********************----------------------------#

def euler_backward(f,y0,x0,xstop,N,args=None):
    h=(xstop-x0)/N
    result=[y0]
    xpar=[x0]
    iter=0
    while iter<N:
        iter+=1
        list_of_args=[f,result[-1],x0,h,args]
        sol=newton_raphson(func,result[-1],0.0001,list_of_args)
        d=call_func(f, args, sol, x0 + h)
        d1=h*d
        d2=d1+result[-1]
        result.append(d2)
        x0+=h
        xpar.append(x0)
    return xpar,result

#-----------------------**********************------------------------#

def rk2(f,y0,x0,xstop,N,args=None):
    h = (xstop - x0) / N
    result=[y0];c=[x0]
    iter=0
    while iter<N:
        iter+=1
        d=call_func(f,args,result[-1],x0)
        k1_by_2=h*d/2
        yn=result[-1]+k1_by_2
        d=call_func(f,args,yn,x0+h/2)
        k2=h*d
        yn=result[-1]+k2
        result.append(yn)
        x0+=h
        c.append(x0)
    return c,result

#------------------------***********************----------------------#

def predictor_corrector(f,y0,x0,xstop,N,args=None):
    h = (xstop - x0) / N
    result = [y0];c=[x0]
    iter=0
    while iter<N:
        iter+=1
        d = call_func(f, args, result[-1], x0)
        k1 = h * d
        yp = result[-1]+ k1 #predicted
        d = call_func(f, args, yp, x0+h)
        k2_p=h* d #predicted k2
        k=(k2_p+ k1)/2 # h(f(y0,x0) +f(y1_p,x0+h)
        yc=result[-1]+ k #corrected
        result.append(yc)
        x0+=h
        c.append(x0)
    return c,result
#--------------------------************************----------------------------#

def rk4(f,y0,x0,xstop,N,args=None):
    h=(xstop-x0)/N
    if type(y0) != list and type(y0) != tuple: y0 = [y0]
    prev_result = y0
    y=[[prev_result[i]] for i in range(len(y0))]
    c = [x0]
    iter=0
    while iter<N:
        iter+=1
        #finding k1/2 and k1/6
        d = call_func(f, args, prev_result, x0)
        if type(d) != list and type(d) != tuple: d = [d]
        k1_by_2 = mul(h / 2, d)
        k1_by_6=mul(h / 6, d)
        #finding k2/2 and k2/3
        yn = add(prev_result, k1_by_2)
        d = call_func(f, args, yn, x0 + h / 2)
        if type(d) != list: d = [d]
        k2_by_2=mul(h/2,d)
        k2_by_3=mul(h/3,d)
        #finding k3
        yn = add(prev_result, k2_by_2)
        d = call_func(f, args, yn, x0 + h / 2)
        if type(d) != list: d = [d]
        k3 = mul(h, d)
        k3_by_3= mul(h/3, d)
        #finding k4
        yn = add(prev_result, k3)
        d = call_func(f, args, yn, x0+h)
        if type(d) != list: d = [d]
        k4_by_6=mul(h/6,d)
        #finding the correct yn
        yn=add(prev_result, k1_by_6, k2_by_3, k3_by_3, k4_by_6)
        prev_result=yn
        for i in range(len(y0)): y[i].append(prev_result[i])
        x0 += h
        c.append(x0)
    return c,y

#----------------------********************-----------------------#



def _shooter_helper(f,ya,yb,a,b,r,s,epsilon,N,args,no_iteration):
    if no_iteration>50: #recursion limit
        print("The limit has reached.")
        return r
    no_iteration+=1
    c1=r[0][-1] #boundry value at choice1
    choice1=r[1][0]
    if abs(c1-yb)<epsilon: return r
    c2=s[0][-1] #boundary value at choice 2
    choice2=s[1][0]
    if abs(c2-yb)<epsilon: return s
    if (c1-yb)*(c2-yb)>0 or c1==c2:
        print("The choices are not good")
        return r
    choice=(yb-c1)*(choice2-choice1)/(c2-c1) + choice1
    q,p=rk4(f,[ya,choice],a,b,N,args)
    c=p[0][-1]
    if abs(c - yb) < epsilon: return p
    if (c1 - yb) * (c - yb) > 0 : #change the choice1
        return _shooter_helper(f, ya, yb, a, b, p, s,epsilon, N, args, no_iteration + 1)
    else: #change the choice2
        return _shooter_helper(f, ya, yb, a, b, r, s,epsilon, N, args, no_iteration + 1)


#shooter mathod
def shooter(f,ya,yb,choice1,choice2,a,b,N,epsilon=0.0001,args=None):
    x,r=rk4(f,[ya,choice1],a,b,N,args)
    x,s=rk4(f,[ya,choice2],a,b,N,args)
    return x, _shooter_helper(f, ya, yb, a, b, r,s ,epsilon, N, args, 0)


#-----------------------**********************-----------------------------#

