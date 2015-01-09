# -*- coding: utf-8 -*-
#/usr/bin/env python
from __future__ import division
from pylab import *

# Parameters for Lennnard-Jones Force
# Contains -> Name:  [characteristic molecular size (nm), binding energy (scaled)]
# Sourced from http://gozips.uakron.edu/~mattice/ps674/lj.html
particles = {'Hydrogen': [.281, 8.6e-3],
             'Flourine': [.283, 52.8e-3],
             'Chlorine': [.335, 173.5e-3]}

def zeroLJForce(sigma):
    """
    Returns the radius where the Lennard-Jones force is zero
    sigma - Characteristic molecular size
    """
    return 2**(1/6) * sigma

def LJForce(energy, sigma, r):
    return 24 * energy * (2 * (sigma**12/r**13) - (sigma**6/r**7))

def magnitude(p):
    """
    Calculate magnitude of the vector
    """
    return sqrt((p**2).sum())

def compute_acceleration(energy, sigma, r):
    """
    Compute the acceleration at aseperation
    energy - Characteristic binding energy
    sigma - Characteristic molecular radius
    r - Seperation vector of Particles
    Treats both particles as identical masses
    """
    return 0.5 * -LJForce(energy, sigma, magnitude(r)) * (r/magnitude(r))

def Leapfrog2D(p1, v1, p2, v2, dt, energy, sigma):
    """
    Leapfrog method in 2D
    p - position vector
    v - velocity vector
    dt - timestep
    energy - Characteristic binding energy
    sigma - Characteristic molecular radius
    returns p, v after dt
    """
    a = compute_acceleration(energy, sigma, p2-p1)
    pn1 = p1 + v1*dt + 0.5 * a * dt**2
    pn2 = p2 + v2*dt - 0.5 * a * dt**2
    a2 = compute_acceleration(energy, sigma, pn2-pn1)
    vn1 = v1 + 0.5 * (a+a2) * dt
    vn2 = v2 - 0.5 * (a+a2) * dt
    return pn1, vn1, pn2, vn2

def EvolveSystem(initial_conditions, tMax, n, energy, sigma):
    """
    Solve the system using the given variables
    initial_condifions - array containing initial positition and velocity as [[p0, v0], [p0,v0]]
    tMax - Time to perform simulation to
    n - number of steps to use
    energy - Characteristic binding energy
    sigma - Characteristic molecular radius
    """
    p1 = zeros([n, 2])
    v1 = zeros_like(p1)
    p2 = zeros([n, 2])
    v2 = zeros_like(p2)
    dt = tMax/n
    p1[0], v1[0], p2[0], v2[0] = initial_conditions
    for i in range(n-1):
        p1[i+1], v1[i+1], p2[i+1], v2[i+1] = Leapfrog2D(p1[i], v1[i], p2[i], v2[i], dt, energy, sigma)
    return p1, v1, p2, v2

if __name__ == '__main__':
    # Plot graphs of Lennard-Jones force for all particles in list
    for molecule in particles:
        sigma, energy = particles[molecule]
        r = linspace(sigma, 1, num=200)
        f = LJForce(energy, sigma, r)
        plot(r, f, label=molecule)

    # Prettify plot
    ylabel('Force')
    xlabel('Separation / nm')
    legend(loc=4)
    title('Lennard-Jones Force for Various Elements')
    ylim(-1.5, 0.5)
    grid()

    # Now do a simulation of two hydroge particles
    figure()
    xlabel('X')
    ylabel('Y')
    grid()
    title('Simulation of two Particles')
    energy, sigma = particles['Hydrogen']
    # Set seperation at 0 LJ force
    initial_position = array([zeroLJForce(sigma), 0])
    # Set veolocity in random direction for KE = binding energy/2
    initial_speed = sqrt(energy)
    angle = 2*pi*random()
    initial_velocity = array([sin(angle)*initial_speed, cos(angle)*initial_speed])
    initial_conditions = vstack((array([0,0]), array([0,0]), initial_position, initial_velocity))

    p1, v1, p2, v2 = EvolveSystem(initial_conditions, 2, 2000, energy, sigma)

    plot(p1[:,0], p1[:,1])
    plot(p2[:,0], p2[:,1])

    # Check Energies
    total_initial_energy = energy/2

    e1 = array([0.5 * magnitude(v)**2 for v in v1])
    e2 = array([0.5 * magnitude(v)**2 for v in v2])
    te = e1+e2

    print('Initial Energy: ', total_initial_energy)
    print('Maximum KE during simulation: ', te.max())
    print('Clearly energy is conserved')

    show()
