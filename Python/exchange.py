import cmath, numpy as np, math, random, compute


def step(r, ext_points, angles, new_ext_point, new_angle):
    d = np.linalg.solve(compute.A(ext_points, angles), compute.a(new_ext_point, new_angle))
    delta = np.array([quot(r[k],d[k]) for k in range(len(d))])
    min_index = np.argmin(delta)
    r = r-delta[min_index]*d
    
    ext_points[min_index] = new_ext_point
    angles[min_index] = new_angle
    r[min_index] = delta[min_index]
    return [r, ext_points, angles]

def quot(x,y):
    if y>0:
        return x/y
    else:
        return math.inf
