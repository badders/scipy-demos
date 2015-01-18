from __future__ import division

from pylab import *
from matplotlib import animation

class NeedlePotential():
    """ Class for animation of potential field.
    Inherit from and override initialPotential and fixPotential at a minimum for
    different initial conditions """
    def __init__(self, N=60, steps=400, iterations_per_frame=4):
        self.N = N
        self.steps = steps
        self.iters_per_frame = iterations_per_frame

        # Create and initialise matplotlib objects
        self.fig = figure()
        self.u = zeros((N,N), dtype=float)
        self.__initialPotential()

    def __initialPotential(self):
        """ Create initial potential grid """

        # Create a 1d array to represent needle
        self.needle = concatenate([zeros(self.N-(self.N//1.2)), ones(self.N//1.2)])
        N = self.N

        # Create the 2d potential
        self.u[-1] = ones(self.N)
        self.u[-1][0] = 0
        self.u[-1][-1] = 0
        self.u[:,N//2] = self.needle

        for i in range(N//4, N - N//4):
            self.u[N//4][i] = 1

        for i in range(N//3, N - N//3):
            self.u[N - N//3][i] = 1

        self.plt = imshow(self.u, origin='lower')
        title('Laplace\'s Equation Solver')
        xticks([])
        yticks([])
        colorbar()

    def __fixPotential(self):
        """ Force potential back to fixed values after clobbering by solver
        """
        N = self.N
        self.u[:,N//2] = self.needle
        for i in range(N//4, N - N//4):
            self.u[N//4][i] = 1

        for i in range(N//3, N - N//3):
            self.u[N - N//3][i] = 1

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
        for x in range(self.iters_per_frame):
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
        self.animator =  animation.FuncAnimation(self.fig, self.animate, repeat=False, frames=self.steps,
                                                 interval=1000/framerate)
        return self.animator

    def save(self, filename):
        self.animator.save(filename)


if __name__ == '__main__':
    solver = NeedlePotential()
    solver.createAnimator()
    #solver.save('laplace.mp4')
    show()
