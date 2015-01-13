from __future__ import division

from pylab import *
from matplotlib import animation

# Conditions for waves
tmax = 60*pi
tstep = tmax/200
xmax = 100
xstep = xmax/300
n = 100 # Number of waves to superimpose
krange = [0.7, 1.3]
wrange = [0.8, 1.2]

# Setup arrays
x = arange(0,xmax,xstep)
t = arange(0,tmax,tstep)

kvalues = linspace(*krange, num=n)
wvalues = linspace(*wrange, num=n)

# Create Matplotlib figure and plot objects
fig = figure()
graphR, = plot([], [], label='Real', color='purple')
graphI, = plot([], [], label='Imaginary', color='orange')
graphA, = plot([], [], label='Absolute', color='black')

# Prettify graph and set information
xlabel('Distance')
ylabel('Signal')
ylim(-1.3,1.3)
xlim(0, xmax)
legend()
grid(True)
title('Wave Packet')

def f(x, t, k=1, w=1, A=1/n):
    """ Wave Function """
    return A * exp(1j * (k*x - w*t))

def init():
    """ Clear graph data """
    graphR.set_data([],[])
    graphI.set_data([],[])
    graphA.set_data([],[])
    return graphR, graphI, graphA

def animate(i):
    """ Adjust paramaters for a frame of animation """
    s = zeros_like(x, dtype=complex)
    for j in range(0, n):
        s += f(x, t[i], kvalues[j], wvalues[j])
    graphR.set_data(x, real(s))
    graphI.set_data(x, imag(s))
    graphA.set_data(x, abs(s))
    return graphR, graphI, graphA

anim = animation.FuncAnimation(fig, animate, init_func=init, frames=int(tmax/tstep), interval=1000/30)

show()
