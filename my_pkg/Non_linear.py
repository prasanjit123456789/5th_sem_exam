def derivative(f, x, args=None, h=0.0001): #derivative
    if args != None: return (f(x + h, args) - f(x - h, args)) / (2 * h)
    return (f(x + h) - f(x - h)) / (2 * h)
#----------------************************--------------------#
def d_derivative(f, x, args=None, h=0.0001): #double derivative
    if args != None: return (f(x + h, args) + f(x - h, args) - 2 * f(x, args))/h ** 2
    return (f(x + h) + f(x - h) - 2 * f(x)) / h ** 2
#----------------************************--------------------#
def polynomial(x,list):
    n = len(list)
    sum = 0
    for i in range(n): sum += list[i] * x ** (n - 1 - i)
    return sum
#-------------*********************************----------------#
def poly_derivative(n,llst,x=0):
    if n==0: return polynomial(x,llst),llst
    else:
        new_list=[]
        l=len(llst)
        for i in range(l-1):
            new_list.append(llst[i]*(l-i-1))
        return poly_derivative(n-1,new_list,x)
#----------------************************--------------------#
def call_func(f,args,*x): #helps in passing the arguments
    if args==None:
        #print("Is not it")
        return f(*x)
    else: return f(*x,args)
#-----------------*******************----------------------#
def make_bracket(f,a,b,factor,args=None):#bracketting the root
    i=0
    if a>b: a,b=b,a
    while i<15 and call_func(f,args,a)*call_func(f,args,b)>0:
        i+=1
        if abs(call_func(f,args,a)) < abs(call_func(f,args,b)): a -= (b - a) * factor
        else: b += (b - a) * factor
    if call_func(f,args,a)*call_func(f,args,b)>0:
        print("Give proper a & b.")
        return None,None
    elif call_func(f,args,a)*call_func(f,args,b)==0:
        if call_func(f,args,a)==0: return a,a
        else: return b,b
    else:
        return a,b
#----------------************************--------------------#
def bisection(func,a,b,factor,epsi=1e-6,args=None): #bisection method and factor is for bracketting
    a,b=make_bracket(func,a,b,factor,args)
    if a==None and b==None: return None  #bracketing was not possible
    elif a==b: return b  #a or b hits the root
    else:
        c=(a+b)/2;i=0
        while(i<200 and b-a>epsi and call_func(func,args,c)!=0):
            i+=1
            if call_func(func,args,c)*call_func(func,args,a)>0: a=c
            else: b=c
            c=(a+b)/2
        if i==50: print("The limit has reached.")
        return c
#----------------************************--------------------#
def regula_falsi(func,a,b,factor,epsi=1e-6,args=None):#regula falsi method
    a,b=make_bracket(func,a,b,factor,args)
    if a==None and b==None: return None #bracketing was not possible
    elif a==b: return b  #a or b hits the root
    else:
        c= (call_func(func,args,b)*a-b*call_func(func,args,a))/(call_func(func,args,b)-call_func(func,args,a));i=0
        #min(abs(c-a),abs(b-a))<epsi to stop loop
        while (i < 200 and min(abs(b - c), abs(a - c))>epsi and call_func(func,args,c) != 0):
            i+=1
            if call_func(func,args,c)*call_func(func,args,a)>0: a=c
            else: b=c
            c = (call_func(func,args,b) * a - b * call_func(func,args,a)) / (call_func(func,args,b) - call_func(func,args,a))
        if i==200: print("The limit has reached.")
        return c
#----------------************************--------------------#
def newton_raphson(f,x0,epsl=1e-6,args=None):
    if call_func(f,args,x0)==0: return x0
    fp=derivative(f,x0,args)
    if fp==0:
        print("Re-enter proper x_0")
        return
    x1 = x0 - call_func(f,args,x0) / fp
    i=0
    while abs(x1-x0)>epsl and call_func(f,args,x1)!=0 and i<200:
        i+=1
        x0=x1
        fp=derivative(f, x0,args)
        if fp == 0:
            print("Re-enter proper x_0")
            return
        x1 = x0 - call_func(f,args,x0) / fp
    return x1
#----------------************************--------------------#
def laguerre(list,a0,epsl=1e-6):
    t = len(list)
    root=[]
    while len(root)<t-1: #no of roots is same as degree
        n=len(list)-1
        if polynomial(a0,list)==0:root.append(a0) # a0 bangs on a root
        else:
            if len(list)==2: #when the polynomial becomes linear
                root.append(-list[1]/list[0])
                return root
            e, f = poly_derivative(1, list, a0)[0], poly_derivative(2, list, a0)[0]
            G=e/polynomial(a0,list)
            H=G**2-f/polynomial(a0,list)
            ll=((n-1)*(n*H-G**2))**0.5
            if abs(G+ll)>abs(G-ll): a=n/(G+ll)
            elif abs(G+ll)==0:
                print("Re-enter proper x_0.")
                return root
            else: a=n/(G-ll)
            a0-=a
            j=0
            while abs(a)>epsl and polynomial(a0,list)!=0 and j<50: #finding the root using G and H
                j+=1
                e,f=poly_derivative(1,list,a0)[0],poly_derivative(2,list,a0)[0]
                G =  e / polynomial(a0,list)
                H = G ** 2 - f / polynomial(a0,list)
                ll = ((n - 1) * (n * H - G**2)) ** 0.5
                if abs(G + ll) > abs(G - ll):a = n / (G + ll)
                elif abs(G + ll) == 0:
                    print("Re-enter proper x_0.")
                    return 0
                else:a = n / (G - ll)
                a0-=a
            if j==50: print("LIMIT HAS REACHED")
            root.append(a0)
        if len(root) != 0: #synthetic division
            new_list = [list[0]]
            last_root = root[len(root) - 1]
            for i in range(1, len(list) - 1):
                new_list.append(last_root * new_list[-1] + list[i])
            list = new_list
    return root
#----------------************************--------------------#

