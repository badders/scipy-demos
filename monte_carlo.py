# -*- coding: utf-8 -*-
from __future__ import division
from pylab import *
from mpl_toolkits.mplot3d import Axes3D

def coneZ(x, y):
    """
    Return Z coordinate of the cone at for corresponding x,y coordinates
    """
    return sqrt(x**2 + y**2)

def withinSphere(x, y, z):
    """
    Retrun whether coordinate is within sphere
    Sphere has radius of 1, and is offset from origin by z+1
    """
    return sqrt(x**2 + y**2 + (z-1)**2) <= 1

def aboveCone(x, y, z):
    """
    Return whether coordinate is aabove the cone
    """
    return coneZ(x, y) < z

if __name__ == '__main__':
    # Create figure with 3D axes
    fig = figure()
    ax = fig.add_subplot(111, projection='3d')

    # Generate and plot Sphere
    u = linspace(0, 2*pi, 100)
    v = linspace(0, pi, 100)

    x = outer(cos(u), sin(v))
    y = outer(sin(u), sin(v))
    z = outer(ones(size(u)), cos(v))+1

    ax.plot_surface(x, y, z, rstride=5, cstride=5, linewidth=0, alpha=0.2, color='b')
    # Reuse x and y and plot cone
    ax.plot_surface(x, y, coneZ(x, y), rstride=5, cstride=5, linewidth=0, alpha=0.4, color='r')

    # Generate random points and check locations
    N = 1500

    points = random([N, 3])*2-1
    points[:,2] += 1

    inside = []
    outside = []

    for p in points:
        if withinSphere(*p) and aboveCone(*p):
            inside.append(p)
        else:
            outside.append(p)

    inside = array(inside)
    outside = array(outside)

    # Plot the random points
    ax.scatter(outside[:,0], outside[:,1], outside[:,2], color='y', label='Outside')
    ax.scatter(inside[:,0], inside[:,1], inside[:,2], color='g', label='Inside')

    ax.set_title('Monte Carlo Volume Finder')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')

    print('Volume: %.2f ' % (8*len(inside)/N))
    show()
