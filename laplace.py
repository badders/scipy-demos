from __future__ import division

from pylab import *
from matplotlib import animation

class NeedlePotential():
    """ Class for animation of potential field.
    Inherit from and overide initialPotential and fixPotential at a minumum for
    different initial conditions """
    def __init__(self, N=40, steps=400):
        self.N = N
        self.steps = steps
        # Create and initialise matplotlib objects
        self.fig = figure()
        self.u = zeros((N,N), dtype=float)
        self.__initialPotential()

    def __initialPotential(self):
        """ Create initial potential grid """
        # Create the 2D array of potential:
        # ___________________
        #         |
        #         |
        #         |
        #         |
        #
        #

        # Create a 1d array to represent needle
        self.needle = concatenate([zeros(self.N-(self.N//1.2)), ones(self.N//1.2)])
        # Create the 2d potential
        self.u[-1] = ones(self.N)
        self.u[-1][0] = 0
        self.u[-1][-1] = 0
        self.u[:,self.N//2] = self.needle
        self.plt = imshow(self.u, origin='lower')
        title('Laplace\'s Equation Solver')
        colorbar()

    def __fixPotential(self):
        """ Force potential back to fixed values after clobbering by solver """
        self.u[:,self.N//2] = self.needle

    def __valuesolve(self, i, j):
        """ Solve for one value """
        self.u[i, j] = 0.25*(self.u[i-1,j] + self.u[i+1,j] + self.u[i,j-1] + self.u[i,j+1])
        self.__fixPotential()

    def solve(self):
        for i  in range(1, self.N-1):
            for j  in range(1, self.N-1):
                self.__valuesolve(i, j)

    def animate(self, i):
        """ Iterate to new data, and draw one animation frame """
        for i  in range(1, self.N-1):
            for j  in range(1, self.N-1):
                self.__valuesolve(i, j)
        hold()
        self.plt = imshow(self.u, origin='lower')
        hold()
        title('Laplace\'s Equation Solver')
        contour(self.u, 6, colors='k')
        return self.plt

    def createAnimator(self, framerate=30):
        self.animator =  animation.FuncAnimation(self.fig, self.animate, repeat=False, frames=self.steps, interval=1000/framerate)
        return self.animator

    def save(self, filename):
        self.animator.save(filename)


if __name__ == '__main__':
    solver = NeedlePotential()
    solver.createAnimator()
    #solver.save('laplace.mp4')
    show()
