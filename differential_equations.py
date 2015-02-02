from __future__ import division

import numpy as np
import matplotlib.pyplot as plt

def magnitude(p):
    """
    Calculate magnitude of the vector
    """
    return np.sqrt((p**2).sum())

def compute_accelerationNewton(p):
    """
    Compute 2D acceleration at position p using Newtonian Mechanics
    """
    r = magnitude(p)
    return - p / r**3

def compute_accelerationEinstein(p):
    """
    Computer 2D acceleration at position p using approximation of General Relativity
    """
    r = magnitude(p)
    return - (p / r**3) - (768 * p / (25 * r**8))

# Set which acceleration function to use for simulations by default
compute_acceleration = compute_accelerationNewton

def Euler2D(p, v, dt):
    """
    Eulers method in 2D
    p - position vector
    v - velocity vector
    dt - timestep
    returns p, v after dt
    """
    a = compute_acceleration(p)
    pn = p + v*dt
    vn = v + a*dt
    return pn, vn

def EulerCromer2D(p, v, dt):
    """
    Euler-Cromer method in 2D
    p - position vector
    v - velocity vector
    dt - timestep
    returns p, v after dt
    """
    a = compute_acceleration(p)
    vn = v + a*dt
    pn = p + vn*dt
    return pn, vn

def Leapfrog2D(p, v, dt):
    """
    Leapfrog method in 2D
    p - position vector
    v - velocity vector
    dt - timestep
    returns p, v after dt
    """
    a = compute_acceleration(p)
    pn = p + v*dt + 0.5 * a * dt**2
    a2 = compute_acceleration(pn)
    vn = v + 0.5 * (a+a2) * dt
    return pn, vn

def OrbitSolve(initial_conditions, tMax, n, method):
    """
    Solve the system using the given variables
    initial_condifions - array containing initial positition and velocity as [p0, v0]
    tMax - Time to perform simulation to
    n - number of steps to use
    method - function to generate new values with
    """
    p = np.zeros([n, 2])
    v = np.zeros_like(p)
    dt = tMax/n
    p[0], v[0] = initial_conditions
    for i in range(n-1):
        p[i+1], v[i+1] = method(p[i], v[i], dt)
    return p, v

def SolveToZero(initial_conditions, dt, method):
    """
    Evolve the system until 2 particles collide
    initial_conditions - array containing initial position and velocity as [p0, v0]
    dt - time step to use
    method - function to generate new values with
    """
    p, v = initial_conditions
    t = 0
    last_p = p*2
    # Loop until particle separation increases (i.e collision has happened)
    while magnitude(p) < magnitude(last_p):
        last_p = p
        p, v = method(p, v, dt)
        t += dt
    return t

if __name__ == '__main__':
    # Setup initial conditions
    r = 10.0
    w = np.sqrt(1 / r**3)
    v = w*r
    p0 = np.array([r, 0.0])
    v0 = np.array([0.0, v])
    i = np.array([p0, v0])
    n = 1000

    # Set time to  2  orbits
    tmax =  4 * np.pi / w
    # Also set up a time for 10 orbits
    tmax10 = 20 * np.pi / w

    # Check the convergence using Newtonian mechanics
    compute_acceleration = compute_accelerationNewton
    pE, _ = OrbitSolve(i, tmax, n, Euler2D)
    pEC, _ = OrbitSolve(i, tmax, n, EulerCromer2D)
    pLF, _ = OrbitSolve(i, tmax, n, Leapfrog2D)

    rE = np.sqrt(pE[:,0]**2 + pE[:,1]**2)
    rEC = np.sqrt(pEC[:,0]**2 + pEC[:,1]**2)
    rLF = np.sqrt(pLF[:,0]**2 + pLF[:,1]**2)

    plt.figure()
    plt.title('Comparison of Interpolation Methods')
    plt.plot(rE, label='Euler')
    plt.plot(rEC, label='Euler-Cromer')
    plt.plot(rLF, label='Leapfrog')
    plt.xlabel('Time')
    plt.ylabel('Radius of Orbit')
    plt.legend(loc='best')

    # Check the convergence using general relativity approximation
    compute_acceleration = compute_accelerationEinstein
    r = 3.0
    w = np.sqrt(1 / r**3)
    v = w*r
    p0 = np.array([r, 0.0])
    v0 = np.array([0.0, v])
    i = np.array([p0, v0])
    T = 2 * np.pi / w
    ts = np.array([5.0, 10, 15])

    for t in ts:
        p, _ = OrbitSolve(i, t, int(t*100), Leapfrog2D)
        rs = np.sqrt(p[:,0]**2 + p[:,1]**2)
        print('T=%d, minimum radius=%f' % (t, rs.min()))
        # Particle orbit radius gets smaller each iteration using this method

    t = SolveToZero(i, T/100, Leapfrog2D)
    print('Particles Collide at t=', t)
    plt.show()