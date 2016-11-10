"""This module implements the multi-scale Lorenz model.

References:
    [1] Lorenz, E. N. (1996). Predictability: A problem partly solved.
        Proc. Seminar on predictability.
    [2] Lorenz, E. N. and Emanuel K. A. (1998). Optimal Sites for
        Supplementary Weather Observations: Simulation with a Small
        Model. Journal of the Atmospheric Sciences.

Author:             Jeremy Shaw
Institution:        Portland State University
Date Created:       17 November 2015
Last Modified Date: 10 January 2016
"""
import numpy as np

import lz96_fortran as fortran

class LZ96(object):
    """Implements the multi-scale model of Edward Lorenz.

    This class models the multi-scale model by Edward Lorenz as in
    [1]. The model consists of the coupled ODEs

    dX[k] / dt = X[k - 1] * (X[k + 1] - X[k - 2]) - X[k] - (h * c / b)
                        * sum_{j = 1}^J Y[j, k] + F                (1)
    dY[j, k] / dt = c * b * Y[j + 1, k] * (Y[j - 1, k] - Y[j + 2, k]) 
                      - c * Y[j, k] + (h * c / b) * X[k]           (2)

    The state vector x is the collection of all variables X_k and Y_jk
    are ordered X_1, X_2, ..., X_K, Y_11, Y_21, ..., Y_J1, Y_12,
    ........, Y_JK.
    
    Attributes:
        b - Set to 10.0 in [1]
        c - Set to 10.0 in [1]
        F - Set to 8.0 in [2]
        h - set to 1.0 in [1]
        J - Number of small-scale Y variables between two X variables
        K - Number of large-scale X variables
        n - Total number of variables in the model
        dt - Time step for the fourth-order Runge-Kutta method
        x0 - Initial state for model integration

    Functionality:
        model - The right side of the coupled ODE system (1)--(2)
        model_tlm - Tangent linear model of the right side of the ODE
            system
        model_adj - The adjoint of the right side of the ODE system
        forecast - Numerical integration of the ODE system via RK4
            method
        forecast_tlm - Tangent linear model of the RK4 method
        forecast_adj - Adjoint model of the RK4 method
    """
    def __init__(self, b=10.0, c=10.0, F=8.0, h=1.0, J=10, K=40, \
                 dt=0.005):
        """Initializes model parameters.

        Initializes the parameters for the Lorenz 96 model. In [1], he
        sets b = c = 10.0, h = 1.0, K = 36, and J is set to 10, but in
        [2], F is set to 8.0 and K (called J there) equals 40.
        """
        self.b = b
        self.c = c
        self.F = F
        self.h = h
        self.J = J
        self.K = K
        self.n = K * (J + 1)
        self.dt = dt

        self.x0 = np.zeros(self.n)
        self.x0[0 : K] = F
        self.x0[(K - 1) / 2] += 0.008

    def model(self, x):
        """The governing equations dX/dt and dY/dt (1)--(2).

        Argument:
            x - The vectorized input [X; Y].

        Returns:
            Model state evaluated at the input.
        """
        return fortran.lz96_model(x, self.b, self.c, self.F, self.h, \
                                  self.J, self.K)

    def model_tlm(self, x, xd):
        """Tangent linear model of the governing equations (1)--(2).

        The tangent linear model performs the product of the Jacobian
        matrix evaluated at x with the direction vector xd.

        Arguments:
            x - State vector
            xd - Direction vector
            
        Returns:
            Jacobian-vector product.
        """
        return fortran.lz96_model_tlm(x, xd, self.b, self.c, self.F, \
                                      self.h, self.J, self.K)

    def model_adj(self, x, yb):
        """Adjoint model of the governing equations (1)--(2).

        The adjoint model performs the product of the transpose of the
        Jacobian evaluated at x with the direction vector yb.

        Arguments:
            x - Vector of x-variables
            yb - Direction vector

        Returns:
            Jacobian transpose-vector product.
        """
        return fortran.lz96_model_adj(x, yb, self.b, self.c, self.F, \
                                      self.h, self.J, self.K)

    def forecast(self, x, ndt=1):
        """RK4 integration of the ode system (1)--(2).

        Argument:
            x - Initial state of the model trajectory
            ndt - Number of time-steps to integrate

        Returns:
            Model forecast of x.
        """
        y = fortran.lz96_forecast(x, self.b, self.c, self.F, self.h, \
                                 self.J, self.K, self.dt)

        for i in xrange(1, ndt):
            y = fortran.lz96_forecast(y, self.b, self.c, self.F, \
                                      self.h, self.J, self.K, self.dt)

        return y

    def forecast_tlm(self, x, xd, ndt=1):
        """Tangent linear model of RK4 integration of (1)--(2).

        Arguments:
            x - The state at which the TLM is initialized
            xd - Direction vector
            ndt - Number of time-steps to integrate

        Returns:
            yd - The derivative of RK4 integration in the direction of
                xd.
        """
        yd = fortran.lz96_forecast_tlm(x, xd, self.b, self.c, self.F,\
                            self.h, self.J, self.K, self.dt)

        for i in xrange(1, ndt):
            x = fortran.lz96_forecast(x, self.b, self.c, self.F, \
                            self.h, self.J, self.K, self.dt)
            yd = fortran.lz96_forecast_tlm(x, yd, self.b, self.c, \
                            self.F, self.h, self.J, self.K, self.dt)

        return yd

    def forecast_adj(self, x, yb, ndt=1):
        """Adjoint model of RK4 integration of (1)--(2).

        Arguments:
            x - The state at which the adjoint model is formed.
            yb - Direction vector
            ndt - Number of time-steps

        Returns:
            xb - The backward integration of yb.
        """
        xlist = [[]] * ndt
        xlist[0] = x
        xb = yb

        # Forward sweep
        for i in xrange(1, ndt):
            xlist[i] = fortran.lz96_forecast(xlist[i - 1], self.b, \
                    self.c, self.F, self.h, self.J, self.K, self.dt)

        # Backward sweep
        for i in xrange(ndt - 1, -1, -1):
            xb = fortran.lz96_forecast_adj(xlist[i], xb, self.b, \
                self.c, self.F, self.h, self.J, self.K, self.dt)
        
        return xb
