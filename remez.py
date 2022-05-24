import cmath, numpy as np, math, random

""" This program uses the refined Remez algorithm invented by Peter Tang, Jan Modersitzki and Bernd Fischer 
to compute complex Chebyshev Polynomials corresponding to discretized domains in the complex plane. """

y_symmetric = False # True if conjugate(E) = E, and conjugate(W) = W
x_symmetric = False # True if -conjugate(E) = E and -conjugate(W) = W
is_even = False 
E = [] # The discrete set
W = [] # The weight function for the discrete set, w(z) for z in E


def main():
    print("hello") 

# compute computes the Chebyshev polynomial corresponding to the set E which has to be initialized via the set_E function
# The variable degree specifies the degree of the chebyshev polynomial to be computed
# threshold denotes the absolute error between the generated linear functional and the generated polynomial. The difference between norms of the Chebyshev polynomial and the generated polynomial lies within this threshold.
def compute(degree, threshold = 0.00000001, iterations = 1000, extremal_set = [], extremal_angles = [], supress_progress = True, ignore_threshold = False):
    global is_even
    is_even = (degree%2 == 0)
    basis = get_basis(degree)
    nbr_of_optimzation_points = get_nbr_of_optimization_points(degree)
    indices, trial_angles = get_init_set(extremal_set, extremal_angles, nbr_of_optimzation_points)
    angles = init_angles(indices, trial_angles, basis)
    coefficients = []
    percentage_counter = 0
    eps = 100000000000000000000000
    print("INITIALIZING RUN...")
    for j in range(iterations):        
        r = np.linalg.solve(A(indices, angles, basis),np.array([1] + [0 for k in range(len(indices)-1)]).T)
        coefficient_data = np.linalg.solve(A(indices, angles, basis).T,c_vector(indices, angles, degree))
        h = coefficient_data[0]
        coefficients = [coefficient_data[k+1] for k in range(len(coefficient_data)-1)]
        p_approx = np.poly1d([0])
        for k in range(len(coefficients)):
            p_approx += coefficients[k]*basis[k]
        max_error, max_index = compute_maximum(p_approx, degree)
        if h>0:
            eps = (max_error-h)/h
        if (j/iterations)*100>percentage_counter*10:
            print(str(percentage_counter*10)+" percentage complete")
            percentage_counter+=1
        if supress_progress==False:
            print("RELATIVE ERROR = "+str(eps))
        if eps<threshold and (ignore_threshold == False):
            print("RELATIVE ERROR < "+str(threshold))
            T_n = np.poly1d([1]+[0 for k in range(degree)])-p_approx
            extremal_set = []
            return [T_n, [E[k] for k in indices], angles, eps]
        else:
            new_index = max_index
            new_angle = cmath.phase(f(E[new_index], degree)-p_approx(E[new_index]))
            r, indices, angles = exchange(r, indices, angles, new_index, new_angle, basis)
            
    print("TIMED OUT, RELATIVE ERROR = "+str(eps))
    T_n = np.poly1d([1]+[0 for k in range(degree)])-p_approx
    extremal_set = []
    return [T_n, [E[k] for k in indices], angles, eps]

def exchange(r, indices, angles, new_index, new_angle, basis):
    d = np.linalg.solve(A(indices, angles, basis), a(new_index, new_angle, basis))
    delta = 10000000000000000000000000000000000
    min_index = 0
    for k in range(len(d)):
        if d[k]>0:
            if delta > r[k]/d[k]:
                delta = r[k]/d[k]
                min_index = k
    r = r-delta*d
    indices[min_index] = new_index
    angles[min_index] = new_angle
    r[k] = delta
    return [r, indices, angles]
    
def get_init_set(extremal_set, extremal_angles, nbr_of_optimization_points):
    E_sample = random.choices(E, k = nbr_of_optimization_points+1)
    while len(extremal_set)<nbr_of_optimization_points:
        for z in E_sample:
            if (z in extremal_set) == False:
                extremal_set.append(z)
                break
    I = np.linspace(0, math.pi, 1000)
    I_sample = random.choices(I, k=nbr_of_optimization_points+1)
    while len(extremal_angles)<nbr_of_optimization_points:
        for t in I:
            if (t in extremal_angles) == False and (t+math.pi in extremal_angles) == False:
                extremal_angles.append(t)
                break
    return [[E.index(z) for z in extremal_set], extremal_angles]

def f(z, degree):
    return z**degree

def get_basis(degree):
    basis = []
    for k in range(degree):
        p = np.poly1d([1]+[0 for l in range(k)])
        if x_symmetric and y_symmetric:
            if is_even and k%2 == 0:
                basis += [p]
            elif is_even == False and (k+1)%2 == 0:
                basis += [p]
        elif y_symmetric:
            basis += [p]
        elif x_symmetric:
            if is_even:
                if k%2 == 0:
                    basis += [p]
                else:
                    basis += [1J*p]
            else:
                if (k+1)%2 == 0:
                    basis += [p]
                else: 
                    basis += [1J*p]
        else:
            basis += [p, 1J*p]
    return basis

def set_E(base_set, weight = [], x_symmetry = False, y_symmetry = False):
    global E, W, x_symmetric, y_symmetric
    E = []
    W = []
    x_symmetric = x_symmetry
    y_symmetric = y_symmetry
    base_set = [z for z in base_set]
    if len(weight) == 0:
        weight = [1 for z in base_set]
    if len(weight) != len(base_set):
        print("Error: Weight and Set do not have the same dimension")
        return
    
    if x_symmetric and y_symmetric:
        for k in range(len(base_set)):
            if np.real(base_set[k])>=0 and np.imag(base_set[k])>=0:
                E.append(base_set[k])
                W.append(weight[k])
    elif x_symmetric:
        for k in range(len(base_set)):
            if np.real(base_set[k])>=0:
                E.append(base_set[k])
                W.append(weight[k])
    elif y_symmetric:
        for k in range(len(base_set)):
            if np.imag(base_set[k])>=0:
                E.append(base_set[k])
                W.append(weight[k])
    else:
        E = base_set
        W = weight

def a(index, alpha, basis):
    return np.array([1]+[np.real(cmath.exp(-1J*alpha)*W[index]*p(E[index])) for p in basis])

def A(indices, angles, basis):
    return np.array([a(indices[k], angles[k], basis) for k in range(len(indices))]).T
    
def c(index, alpha, degree):
    return np.real(cmath.exp(-1J*alpha)*W[index]*f(E[index], degree))

def c_vector(indices, angles, degree):
    return np.array([c(indices[k], angles[k], degree) for k in range(len(indices))]).T

def get_nbr_of_optimization_points(degree):
    if x_symmetric and y_symmetric:
        return int(degree/2)+1 
    elif x_symmetric == False and y_symmetric == False:
        return 2*degree + 1
    else:
        return degree + 1

def compute_maximum(p, degree):
    temp_max = 0
    temp_index = 0
    for index in range(len(E)):
        error = np.abs(W[index]*(f(E[index], degree)-p(E[index])))
        if error > temp_max:
            temp_max = error
            temp_index = index
    return [temp_max, temp_index]

def init_angles(indices, trial_angles, basis):
    r = np.linalg.solve(A(indices,trial_angles, basis), np.array([[1]+[0 for k in range(len(indices)-1)]]).T)
    angles = [s for s in trial_angles]
    for k in range(len(r)):
        if r[k]<0:
            angles[k]+=math.pi
    return angles



if __name__=="__main__":
    main()

