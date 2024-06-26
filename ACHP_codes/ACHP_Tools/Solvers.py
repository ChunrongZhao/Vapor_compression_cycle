from __future__                         import division, print_function, absolute_import
import numpy                            as np
from scipy.linalg                       import inv
from numpy                              import array, dot
from scipy.optimize                     import fsolve

# -----------------------------------------------------------------------------------
def MultiDimNewtRaph(f, x0, dx=1e-6, args=(), ytol=1e-5, w=1.0, JustOneStep=False):
    """
    A Newton-Raphson solver where the Jacobian is always re-evaluated rather than
    re-using the information as in the fsolve method of scipy.optimize
    """
    # Convert to numpy array, force type conversion to float
    x                                   = np.array(x0, dtype=float)
    error                               = 999
    J                                   = np.zeros((len(x), len(x)))

    # If a float is passed in for dx, convert to a numpy-like list the same shape as x
    if isinstance(dx, int):
        dx                              = dx * np.ones_like(float(x))
    elif isinstance(dx, float):
        dx                              = dx * np.ones_like(x)

    r0                                  = array(f(x, *args))
    while abs(error) > ytol:
        # Build the Jacobian matrix by columns
        for i in range(len(x)):
            epsilon                     = np.zeros_like(x)
            epsilon[i]                  = dx[i]
            J[:, i]                     = (array(f(x + epsilon, *args)) - r0) / epsilon[i]
        # https://numpy.org/doc/stable/reference/generated/numpy.linalg.inv.html
        # inv: compute the inverse of a matrix
        v                               = np.dot(-inv(J), r0)
        x                               = x + w * v
        # Calculate the residual vector at the new step
        r0                              = f(x, *args)
        error                           = np.max(np.abs(r0))
        # Just do one step and stop
        if JustOneStep == True:
            return x
    return x


def Broyden(f, x0, dx=1e-5, args=(), ytol=1e-5, w=1.0, itermax=10, JustOneStep=False):
    """
    Broyden's method:
    Generalization of the secant method to nonlinear systems.
    Roughly speaking, the secant method replaces the derivative
    by a finite difference approximation
    """
    x0                                  = np.array(x0, dtype=float)
    x1                                  = x0.copy()
    error                               = 999
    A1                                  = np.zeros((len(x0), len(x0)))
    A0                                  = np.zeros((len(x0), len(x0)))

    # If a float is passed in for dx, convert to a numpy-like list the same shape as x
    if isinstance(dx, float):
        dx                              = dx * np.ones_like(x0)

    F0                                  = array(f(x0, *args))
    iter                                = 1
    x1                                  = x0.copy()
    F1                                  = F0.copy()
    while abs(error) > ytol:
        if iter == 1:
            # Build the Jacobian matrix by columns
            for i in range(len(x0)):
                epsilon                 = np.zeros_like(x0)
                epsilon[i]              = dx[i]
                A0[:, i]                = (array(f(x0 + epsilon, *args)) - F0) / epsilon[i]
            # Get the difference vector
            x1                          = x0 - np.dot(inv(A0), F0)
            # Just do one step and stop
            if JustOneStep == True:
                return x1
            iter                        += 1

        elif iter > 1 and iter < itermax:
            # Jacobian updating parameters
            S                           = x1 - x0
            d                           = np.dot(S.T, S)
            F1                          = array(f(x1, *args))
            Y                           = F1 - F0
            A1                          = A0 + 1.0 / d * dot((Y - dot(A0, S)), S.T)
            x2                          = x1 - np.dot(inv(A1), F1)
            # Update values
            x0                          = x1
            x1                          = x2
            F0                          = F1
            A0                          = A1
            error                       = np.sqrt(np.sum(np.power(F1, 2)))
            iter                        += 1
        else:
            return np.nan*np.ones_like(x0)
    return x1


def DE(fobj, bounds, mut=0.8, crossp=0.7, popsize=20, its=1000):
    """
    Differential Evolutionary (DE) Algorithm:
    R. Storn and K. Price "Differential Evolution - A Simple and Efficient Heuristic
    for global Optimization over Continuous Spaces", Journal of Global Optimization, 1997
    https://pablormier.github.io/2017/09/05/a-tutorial-on-differential-evolution-with-python/
    """
    dimensions                          = len(bounds)
    pop                                 = np.random.rand(popsize, dimensions)
    min_b, max_b                        = np.asarray(bounds).T
    diff                                = np.fabs(min_b - max_b)
    pop_denorm                          = min_b + pop * diff
    fitness                             = np.asarray([fobj(ind) for ind in pop_denorm])
    best_idx                            = np.argmin(fitness)
    best                                = pop_denorm[best_idx]

    for i in range(its):
        for j in range(popsize):
            idxs                        = [idx for idx in range(popsize) if idx != j]
            a, b, c                     = pop[np.random.choice(idxs, 3, replace=False)]
            # https://numpy.org/doc/stable/reference/generated/numpy.clip.html
            # Clip (limit) the values in an array.
            mutant                      = np.clip(a + mut * (b - c), 0, 1)
            cross_points                = np.random.rand(dimensions) < crossp

            if not np.any(cross_points):
                cross_points[np.random.randint(0, dimensions)] = True
            trial                       = np.where(cross_points, mutant, pop[j])
            trial_denorm                = min_b + trial * diff
            f                           = fobj(trial_denorm)
            if f < fitness[j]:
                fitness[j]              = f
                pop[j]                  = trial
                if f < fitness[best_idx]:
                    best_idx            = j
                    best                = trial_denorm
        # https://blog.csdn.net/mieleizhi0522/article/details/82142856
        yield best, fitness[best_idx]


# -----------------------------------------------------------------------------------
if __name__=='__main__':

    def f(x):
        print('.', end='')
        return [1 - 4 * x[0] + 2 * x[0]**2 - 2 * x[1]**3, -4 + x[0]**4 + 4 * x[1] + 4 * x[1]**4]

    _x          = Broyden(f, [0.1, 0.7], ytol=1e-8)
    print(_x)
    print(f(_x))

    _x          = fsolve(f, [0.1, 1.0])
    print(_x)
    print(f(_x))
