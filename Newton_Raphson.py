#!/usr/bin/env python3

from scipy import optimize
import numpy as np


def f(x):
    return (x**3 - 1)


def iterative_improve(good_enough, improve):

    """ """
    def iterate(guess):
        next_try = improve(guess)
        if good_enough(guess, next_try):
            return next_try
        else:
            return iterate(next_try)
    return lambda x: iterate(x)


def fixed_point(f, first_guess, epsilon=1.0e-12):
    # def good_enough(v1, v2):
    #    return abs(v1-v2)/v2 < epsilon
    return iterative_improve(good_enough, f)(first_guess)


def good_enough(v1, v2, epsilon=1e-8):
    return abs((v1 - v2) / v2) < epsilon


def approx_deriv(f, dx=1.0e-8):
    return lambda x: float(f(x + dx) - f(x)) / dx


def newton_transform(f):
    return lambda x: x - float(f(x)) / (approx_deriv(f)(x))


def newton_raphson(f, guess):
    """
    newton_raphson(f, guess)

    Newton Raphson procedure

    Parameters
    ----------
    f : function

    Returns
    -------
    float
        Zero of the specified function.
    """
    return fixed_point(newton_transform(f), guess)


def fNikuradse(fw, Rew):
    """Nikuradse equation modified by for solving in a Newton-Raphson scheme.

    """
    nikuradse = (fw * (0.86 * np.log(4 * Rew * np.sqrt(fw)) - 0.8) ** 2
                 - 1)
    return nikuradse


def fChezy_wall(xRef):
    """Computes the Chezy friction coefficient for the wall region"""
    # Do a first estimate of fw0 with the Blasius equation, as modified by
    # Chiew and Parker in 1994
    fw0 = 0.301 * xRef ** 0.2
    Rew0 = fw0 / xRef
    convergence = False
    while not convergence:
        Rew1 = Rew0
        fw0 = newton_raphson(lambda fw: fNikuradse(fw, Rew0), 0.01)
        Rew0 = fw0 / xRef
        convergence = good_enough(Rew1, Rew0)
    fw = fw0 / 8
    return fw, Rew0

# fprime is not provided so we use Secant method
root = optimize.newton(f, 1.5)

# fprime2 is provided so Halley's method is used.
root1 = optimize.newton(f, 1.5, fprime2=lambda x: 6*x)

# Only fprime is provided so Newton-Raphson is used:
root2 = optimize.newton(f, 1.5, fprime=lambda x: 3 * x**2)

# fprime1 and fprime2 are provided so Halley's method is used
root3 = optimize.newton(f, 1.5, fprime=lambda x: 3 * x**2,
                        fprime2=lambda x: 6*x)

root4 = newton_raphson(f, 1.5)


def main():
    print(root, root1, root2, root3, root4)
    return


if __name__ == '__main__':
    main()
