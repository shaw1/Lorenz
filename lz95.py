"""This module implements the Lorenz 40-variable model.

References:
    [1] Lorenz, E. N. (1996). Predictability: A problem partly solved. Proc.
        Seminar on predictability.
    [2] Lorenz, E. N. and Emanuel K. A. (1998). Optimal Sites for
        Supplementary Weather Observations: Simulation with a Small Model.
        Journal of the Atmospheric Sciences.

Author:             Jeremy Shaw
Institution:        Portland State University
Date Created:       17 November 2015
Last Modified Date: 20 November 2015
"""

import lz95_fortran as fortran

class LZ95(object):
    """Implements the model of Edward Lorenz in [2].

    This class implements the the model given by the coupled differential
    equations dX[k] / dt = X[k] * (X[k + 1] - X[k - 2]) - X[k] + F.

    Attributes:
        F - Set to 8.0 in [2]
        n - Number of states in the model: set to 40 in [2]
        dt - Time step for the fourth-order Runge-Kutta method
        x0 - Initial state for model integration

    Functionality:
        model - The right side of the coupled ODE system
        model_tlm - Tangent linear model of the right side of the ODE system
        model_adj - The adjoint of the right side of the ODE system
        forecast - Numerical integration of the ODE system via RK4 method
        forecast_tlm - Tangent linear model of the RK4 method
        forecast_adj - Adjoint model of the RK4 method
    """
    def __init__(self, F=8.0, n=40, dt=0.05):
        """Initializes model parameters.

        Initializes the parameters for the Lorenz 95 model. In [2], he sets
        F is set to 8.0 and n = 40. The initial state x0 is set as in [2] as
        well.
        """
        from numpy import zeros
        
        self.F = F
        self.n = n
        self.dt = dt

        self.x0 = zeros(n)
        self.x0[0 : n] = F
        self.x0[(n - 1) / 2] += 0.008

    def model(self, x):
        """Models the right side of the governing equations dX/dt.

        Argument:
            x - The vectorized input

        Returns:
            Model state evaluated at the input.
        """
        return fortran.lz95_model(x, self.F)

    def model_tlm(self, x, xd):
        """Tangent linear model of the govering equations dX/dt.

        The tangent linear model performs the product of the Jacobian matrix
        evaluated at x with the direction vector xd.        
        
        Arguments:
            x - Vector of x-variables
            xd - Direction vector

        Returns:
            Jacobian-vector product.
        """
        return fortran.lz95_model_tlm(x, xd, self.F)

    def model_adj(self, x, yb):
        """Adjoint model of the governing equations dX/dt.

        The adjoint model performs the product of the transpose of the
        Jacobian evaluated at x with the direction vector yb.

        Arguments:
            x - Vector of x-variables
            yb - Direction vector

        Returns:
            Jacobian transpose-vector product.
        """
        return fortran.lz95_model_adj(x, yb, self.F)

    def forecast(self, x, ndt=1):
        """RK4 integration of the ode system dX/dt.

        Argument:
            x - Initial state of the model trajectory
            ndt - Number of time-steps to integrate

        Returns:
            Model forecast of x.
        """
        y = fortran.lz95_forecast(x, self.F, self.dt)

        for i in xrange(1, ndt):
            y = fortran.lz95_forecast(y, self.F, self.dt)

        return y

    def forecast_tlm(self, x, xd, ndt=1):
        """Tangent linear model of RK4 integration.

        Arguments:
            x - The state at which the TLM is initialized
            xd - Direction vector
            ndt - Number of time-steps to integrate

        Returns:
            yd - The derivative of RK4 integration in the direction of xd.
        """
        yd = fortran.lz95_forecast_tlm(x, xd, self.F, self.dt)

        for i in xrange(1, ndt):
            x = fortran.lz95_forecast(x, self.F, self.dt)
            yd = fortran.lz95_forecast_tlm(x, yd, self.F, self.dt)

        return yd

    def forecast_adj(self, x, yb, ndt=1):
        """Adjoint model of RK4 integration.

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
            xlist[i] = fortran.lz95_forecast(xlist[i - 1], self.F, self.dt)

        # Backward sweep
        for i in xrange(ndt - 1, -1, -1):
            xb = fortran.lz95_forecast_adj(xlist[i], xb, self.F, self.dt)
        
        return xb
