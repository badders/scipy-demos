from __future__ import division

from pylab import *
from matplotlib import animation

def f(x, t, k=1, w=1, A=1):
    return A * e ** (1j * (k*x - w*t))

tmax = 40*pi
tstep = 0.5

x = arange(0,200,0.1)
t = arange(0,tmax,tstep)

k1 = 1
w1 = 1
k2 = 1.1
w2 = 1.2

fig = figure()
graphR, = plot([], [], label='Real')
graphI, = plot([], [], label='Imaginary', color='r')
graphA, = plot([], [], label='Absolute', color='k')

xlabel('Distance')
ylabel('Signal')
ylim(-2.5,2.5)
xlim(0, 200)
legend()
grid(True)
title('Funky Beats')

def init():
    graphR.set_data([],[])
    graphI.set_data([],[])
    graphA.set_data([],[])
    return graphR, graphI, graphA

def animate(i):
    s = f(x, t[i], k=k1, w=w1)+f(x, t[i], k=k2, w=w2)
    graphR.set_data(x, real(s))
    graphI.set_data(x, imag(s))
    graphA.set_data(x, abs(s))
    return graphR, graphI, graphA

anim = animation.FuncAnimation(fig, animate, init_func=init, frames=int(tmax/tstep), interval=1000/30)

show()

"""
Keeping w the same for both waves, the beats do not move. Keeping k the same for both waves, you dont get
beats, just a constant amplitude at all distance oscillating with time.

When the velocites differ you get funky wave patterns under the absolute wave
"""
