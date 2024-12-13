import cmath, numpy as np, math, random, data, compute, exchange
import matplotlib.pyplot as plt

def run(tolerance = 1e-8, iterations = 60, show_step=False,show_error = False):
    
    nbr_optimisation_pts = len(data.basis(1))+1
    extremal_set, trial_angles = get_init_set(nbr_optimisation_pts)
    angles = init_angles(extremal_set, trial_angles)
    coefficients = []
    errors = []
    for j in range(iterations):        
        functional_weight = np.linalg.solve(compute.A(extremal_set, angles),np.array([1] + [0 for k in range(nbr_optimisation_pts-1)]).T)
        functional_data = np.linalg.solve(compute.A(extremal_set, angles).T,compute.c_vector(extremal_set, angles))
        functional_value = functional_data[0]
        coefficients = [functional_data[k+1] for k in range(len(functional_data)-1)]
        
        deviation, ext_point = compute.maximum(coefficients)
        if(show_step):
            plt.plot(np.linspace(0,1,1000),[-compute.neg_absolute_difference(x,coefficients) for x in np.linspace(0,1,1000)])
            plt.plot([ext_point],[-compute.neg_absolute_difference(ext_point,coefficients)],'o')
            plt.show()
        if functional_value>0:
            eps = (deviation-functional_value)/functional_value
            errors.append(eps)
            print('Relative error = ' + str(eps))
            if eps<tolerance:
                print('success!')
                if(show_error):
                    plt.yscale('log')
                    plt.plot(range(len(errors)),errors)
                    plt.show()
                return [coefficients, extremal_set, angles, eps]
                
        
        ext_angle = compute.phase(ext_point, coefficients)
        functional_weight, extremal_set, angles = exchange.step(functional_weight, extremal_set, angles, ext_point, ext_angle )
    print('algorithm time-out')
    if(show_error):
        plt.yscale('log')
        plt.plot(range(len(errors)),errors)
        plt.show()
    return [coefficients, extremal_set, angles, eps]

def get_init_set(nbr_of_optimization_points):
    cheb_nodes = [np.cos(((2*(k+1)-1)*math.pi)/(2*nbr_of_optimization_points))/2+1/2 for k in range(nbr_of_optimization_points)]
    angles = random.sample(np.linspace(0,math.pi,10000).tolist(), nbr_of_optimization_points)
    return [cheb_nodes, angles]
    
def init_angles(ext_points, trial_angles):
    functional_weight = np.linalg.solve(compute.A(ext_points,trial_angles), np.array([[1]+[0 for k in range(len(ext_points)-1)]]).T)
    angles = [s for s in trial_angles]
    for k in range(len(functional_weight)):
        if functional_weight[k]<0:
            angles[k]+=math.pi
    return angles


