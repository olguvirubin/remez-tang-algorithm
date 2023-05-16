import cmath, numpy as np, math, random
#This is the function that should be altered by the user
#The function f is the function to approximate, this should be f composed with gamma which is parametrised from 0 to 1
#basis should give the value of the basis functions in vector form
#w should return the value of the weight from 0 to 1
""" def f(z):
    return z**(10)

def basis(z):
    return [z**k for k in range(10)]

def w(z):
    return 1 """

#Example Bernoulli Lemniscate of odd degree
degree = 61
def complete_gamma(t):
    if 0<t<=1/4:
        return gamma(t*4)
    elif 1/4<t<=1/2:
        return -np.real(gamma(2-4*t))+1J*np.imag(gamma(2-4*t))
    elif 1/2<t<=3/4:
        return -np.real(gamma(4*(t-1/2)))-1J*np.imag(gamma(4*(t-1/2)))
    elif 3/4<t<=1:
        return np.real(gamma(4-4*t))-1J*np.imag(gamma(4-4*t))

def gamma(t):
    return np.sqrt(2*np.cos(t*math.pi/2))*cmath.exp(1J*t*math.pi/4)

def f(t):
    return gamma(t)**degree

def basis(t):
    return [gamma(t)**(2*k+1) for k in range(int((degree-1)/2))]

def w(t):
    return 1

def get_chebyshev(coefficients):
    P = np.poly1d([1]+[0 for k in range(int(degree))])
    for k in range(len(coefficients)):
        P = P+np.poly1d([-coefficients[k]]+[0 for k in range(2*k+1)])
    return P

""" #Example weight (1+z)^s on unit disk

def gamma(t):
    return cmath.exp(1J*t*math.pi)

def f(t):
    return gamma(t)**21

def basis(t):
     return [gamma(t)**k for k in range(21)]

def w(t):
    s = 1/2
    return np.sqrt(np.power(np.abs(1+gamma(t)),s))

def get_chebyshev(coefficients, tolerance):
    P = np.poly1d([1]+[0 for k in range(21)])
    for k in range(len(coefficients)):
        P = P+np.poly1d([-coefficient_tolerance_adjusted(coefficients[k], tolerance)]+[0 for k in range(k)])
    return P

def coefficient_tolerance_adjusted(coefficient, tolerance):
    if np.abs(coefficient)<10*tolerance:
        return 0
    else: 
        return coefficient

 """

""" #Example m-star
m = 5
n = 11
l = 3
def gamma(t):
    return  t*np.power(2,1/m)

def complete_gamma(t):
    return cmath.exp(2*1J*math.pi*np.floor(t*2*m)/(2*m))*2*m*(t-np.floor(t*2*m)/(2*m))*np.power(2,1/m)

def f(t):
    return (gamma(t))**(2*n*m+l)

def basis(t):
    return [(gamma(t))**(2*k*m+l) for k in range(n)]

def w(t):
    return 1

def get_chebyshev(coefficients):
    P = np.poly1d([1]+[0 for k in range(2*n*m+l)])
    for k in range(len(coefficients)):
        P = P+np.poly1d([-coefficients[k]]+[0 for k in range(2*k*m+l)])
    return P """
 
""" #Example hypo-cycloid
m=9
n=9
l=2

def complete_gamma(t):
    return cmath.exp(1J*2*math.pi*t)+cmath.exp(-1J*2*math.pi*t*(m-1))/(m-1)

def gamma(t):
    return cmath.exp(1J*2*math.pi*t/m)+cmath.exp(-1J*2*math.pi*t*(m-1)/m)/(m-1)

def f(t):
    return (gamma(t))**(m*n+l)

def basis(t):
    return [(gamma(t))**(m*k+l) for k in range(n)]

def w(t):
    return 1

def get_chebyshev(coefficients):
    P = np.poly1d([1]+[0 for k in range(m*n+l)])
    for k in range(len(coefficients)):
        P = P+np.poly1d([-coefficients[k]]+[0 for k in range(m*k+l)])
    return P """