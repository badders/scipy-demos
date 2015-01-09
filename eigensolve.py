from __future__ import division
from pylab import *
from scipy import integrate

# Physical Constants
eV = 1.60217646e-19
hb = 1.05457148e-34
m = 9.10938188e-31

def fequal(f1, f2, sigma):
    """
    Check two floating point numbers are within sigma of each other
    """
    return abs(f1-f2) < sigma

def classicalLimit(n):
    """
    Return classical limits of potential field for energy level n
    Equation from PX2131 lecture notes
    """
    return sqrt((2*n+1)*hb / (m*w))

def V(x):
    """ Create potential well """
    return 0.5 * x**2 * m * w**2

def expectedEigenValue(n, w):
    """" Return theoreticallly expected eigenvalue
    n - Quantum number
    w - Omega
    """
    return ((n+0.5) * hb * w)/eV

def iterateWave(E, dx, Vx):
    """
    Iterate out the wave function for given parameters.
    E - energy (eV)
    dx - step size
    Vx - potential field
    """
    N = len(Vx)
    g = -(dx**2) * 2 * m * (E*eV-Vx) / hb**2
    psi = zeros(N, dtype=float64)
    psi[1] = Vx[0]/1e6
    for i in range(1,N-1):
        psi[i+1]=2*psi[i]-psi[i-1]+g[i]*psi[i]
    return psi

def between(a, b):
    """ Helper function to check 0 is between two values """
    return max(a, b) > 0 and min(a, b) < 0

def normaliseFunction(psi, dx):
    """
    Normalise the eigenfunction. i.e. integral psi**2 = 1
    """
    r = integrate.simps(psi**2, dx=dx)
    factor = sqrt(1 / r)
    return psi * factor

def search(E1, E2, Vx, dx):
    """
    Trinary homing algorithm around E searching for correct eigenvalue
    E1, E2 - Range of energys to search within (eV)
    Vx - potential field
    dx - step size of potential array
    """
    # Dont care about more than 9 sig figs
    if fequal(E1,E2,1e-10):
        return E1
    # Setup search space
    # FIXME: can parallelise this by increasing number of search spaces
    Es = linspace(E1, E2, num=3)
    rs = []
    for Ei in Es:
        psi = iterateWave(Ei, dx, Vx)
        rs.append([Ei, psi[-1]])
    # Set search boundry based on maxima in first 2/3 of psi data
    sigma = psi[:2*N//3].max() / 1e6
    for E, r in rs:
        if fabs(r) <= sigma:
            return E
    # Check which region true value lies in
    if between(rs[0][1], rs[1][1]):
        return search(rs[0][0], rs[1][0], Vx, dx)
    elif between(rs[1][1], rs[2][1]):
        return search(rs[1][0], rs[2][0], Vx, dx)
    else:
        raise ValueError('No eigen value in range %f %f' % (E1, E2))

if __name__ == '__main__':
    # Parameters for program
    w = 1e15 * 2 * pi # omega
    N = 5000 # Resolution of search space

    # Established approximates by hand, used as starting points in search
    estimates = {0: [1, 3], 1: [5, 7], 2: [9, 11], 3:[12,15], 4:[16,20]}

    print('Quantum Number\tEigenvalue (eV)\tExpected Eigenvalue\tAccuracy')
    r = []

    # FIXME: Could run these in parallel, no need to do sequentially
    for n in range(5):
        limit = classicalLimit(n)
        x = linspace(-4*limit, 4*limit, num=N)
        dx = 8 * limit / N
        xs = x / sqrt(hb / (m*w))
        Vx = V(x)
        E = search(estimates[n][0], estimates[n][1], Vx, dx)
        r.append((n, xs, dx, iterateWave(E, dx, Vx)))
        expected = expectedEigenValue(n, w)
        print('%d\t\t%.9f\t%.9f\t\t%.4f%%' % (n, E, expected, 100-100*abs(E-expected)/expected))

    xlimit = abs(r[-1][1][0])
    for n, xs, dx, psi in r:
        # Trim end of psi
        t = n*300+1
        psi = psi[:-t]
        xs = xs[:-t]
        psi = normaliseFunction(psi, dx)
        # Extend each function out to the maximum of all
        xs = concatenate(([-xlimit], xs, [xlimit]))
        psi = concatenate(([0], psi, [0]))
        plot(xs, psi, label='n=%d' % n, linewidth=2)

    xlim(-5.5, 5.5)
    grid(True)
    legend(loc='best')
    title('Eigen Functions')
    xlabel('x / $\sqrt{\hbar / (m\omega)}$')
    ylabel('$\Psi$')
    show()
