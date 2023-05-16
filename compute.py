import cmath, numpy as np, math, random, data, scipy.optimize as sp

def a(ext_point, angle):
    return np.array([1]+[np.real(cmath.exp(-1J*angle)*z*data.w(ext_point)) for z in data.basis(ext_point)])

def A(ext_points, angles):
    return np.array([a(ext_points[k], angles[k]) for k in range(len(ext_points))]).T
    
def c(ext_point, angle):
    return np.real(cmath.exp(-1J*angle)*data.w(ext_point)*data.f(ext_point))

def c_vector(ext_points, angles):
    return np.array([c(ext_points[k], angles[k]) for k in range(len(ext_points))]).T

def maximum(coefficients):
    bnds = [(0,1)]
    x0 = dense_sampling(coefficients)
    res = sp.minimize(neg_absolute_difference, x0=x0, bounds=bnds, args = (coefficients),tol=1e-10)
    return [-res.fun,res.x[0]]

def phase(extremal_point, coefficients):
    return cmath.phase(difference(extremal_point, coefficients))

def difference(x, coefficients):
    return data.w(x)*(data.f(x)-np.dot(np.array(coefficients),np.array(data.basis(x))))

def dense_sampling(coefficients):
    I = np.linspace(0,1,1000)
    max_val = -math.inf
    ext_point = 0
    for x in I:
        temp_val = np.abs(difference(x, coefficients))
        if temp_val>max_val:
            max_val = temp_val 
            ext_point = x
    return ext_point

def neg_absolute_difference(x, coefficients):
    return -np.abs(difference(x, coefficients))