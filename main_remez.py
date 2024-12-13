import cmath, numpy as np, math, random, data, algorithm, compute, scipy.optimize as sp
import matplotlib.pyplot as plt

def main():
    coefficients = algorithm.run(iterations=200, tolerance=1e-15, show_step=False, show_error=False)[0]
    
    I = np.linspace(0,1,10000)
    plt.plot(I, [-compute.neg_absolute_difference(x,coefficients) for x in I])
    plt.show()
    plt.plot(I, [-compute.neg_absolute_difference(x,coefficients) for x in I])
    plt.savefig('absolute_val.png', dpi=600)
    plt.close()
    T = data.get_chebyshev(coefficients)
    print(T)
    
    
    plt.plot([np.real(data.complete_gamma(t)) for t in I],[np.imag(data.complete_gamma(t)) for t in I])
    plt.plot(np.real(T.r), np.imag(T.r),'o')
    plt.gca().set_aspect('equal', adjustable='box')
    plt.show()
    plt.plot([np.real(data.complete_gamma(t)) for t in I],[np.imag(data.complete_gamma(t)) for t in I])
    plt.plot(np.real(T.r), np.imag(T.r),'o')
    plt.gca().set_aspect('equal', adjustable='box')
    plt.savefig('domain.png', dpi=600)
    plt.close()


#new comment yes
if __name__ == "__main__":
    main()