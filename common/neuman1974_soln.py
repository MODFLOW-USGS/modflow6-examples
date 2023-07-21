import numpy as np
from math import sqrt
import matplotlib.pyplot as plt

# Find a root of a function using Brent's method within a bracketed range
from scipy.optimize import brentq

# Solve definite integral using Fortran library QUADPACK
from scipy.integrate import quad

# Zero Order Bessel Function
from scipy.special import j0, jn_zeros

__all__ = ["RadialUnconfinedDrawdown"]

pi = 3.141592653589793
sin = np.sin
cos = np.cos
sinh = np.sinh
cosh = np.cosh
exp = np.exp


def _find_hyperbolic_max_value():
    seterr = np.seterr()
    np.seterr(all="ignore")
    inf = np.inf
    x = 10.0
    delt = 1.0
    for i in range(1000000):
        x += delt
        try:
            if inf == sinh(x):
                break
        except:
            break
    np.seterr(**seterr)
    return x - delt


_hyperbolic_max_value = _find_hyperbolic_max_value()


def _find_hyperbolic_equivalent_value():
    x = 10.0
    delt = 0.0001
    for i in range(1000000):
        x += delt
        if x > _hyperbolic_max_value:
            break
        try:
            if sinh(x) == cosh(x):
                return x
        except:
            break
    return x - delt


_hyperbolic_equivalence = _find_hyperbolic_equivalent_value()


class RadialUnconfinedDrawdown:
    """
    Solves the drawdown that occurs from pumping from partial penetration
    in an unconfined, radial aquifer. Uses the method described in:
    Neuman, S. P. (1974). Effect of partial penetration on flow in
    unconfined aquifers considering delayed gravity response.
    Water resources research, 10(2), 303-312.
    """

    hyperbolic_max_value = _hyperbolic_max_value
    hyperbolic_equivalence = _hyperbolic_equivalence

    bottom: float
    Kr: float
    Kz: float
    Ss: float
    Sy: float
    well_top: float
    well_bot: float
    saturated_thickness: float

    _sigma: float
    _beta: float

    def __init__(
        self,
        bottom_elevation,
        hydraulic_conductivity_radial=None,
        hydraulic_conductivity_vertical=None,
        specific_storage=None,
        specific_yield=None,
        well_screen_elevation_top=None,
        well_screen_elevation_bottom=None,
        water_table_elevation=None,
        saturated_thickness=None,
    ):
        """
        Initialize unconfined, radial groundwater model to solve drawdown
        at an observation location in response to pumping at the center of
        the model (that is, the well extracts water at radius = 0).

        Parameters
        ----------
        rad : int
            radial band number (0 to nradial-1)

        bottom_elevation : float
            Elevation of the impermeable base of the model ($L$)
        hydraulic_conductivity_radial : float
            Radial direction hydraulic conductivity of model ($L/T$)
        hydraulic_conductivity_vertical : float
            Vertical (z) direction hydraulic conductivity of model ($L/T$)
        specific_storage : float
            Specific storage of aquifer ($1/T$)
        specific_yield : float
            Specific yield of aquifer ($-$)
        well_screen_elevation_top : float
            Pumping well's top screen elevation ($L$)
        well_screen_elevation_bottom : float
            Pumping well's bottom screen elevation ($L$)
        water_table_elevation : float
            Initial water table elevation. Note, saturated_thickness (b) is
            calculated as $water_table_elevation - bottom_elevation$ ($L$)
        saturated_thickness : float
            Specify the initial saturated thickness of the unconfined aquifer.
            Value is used to calculate the water_table_elevation. If
            water_table_elevation is defined, then saturated_thickness input
            is ignored and set to
            $water_table_elevation - bottom_elevation$ ($L$)
        """

        self.bottom = float(bottom_elevation)
        self.Kr = self._float_or_none(hydraulic_conductivity_radial)
        self.Kz = self._float_or_none(hydraulic_conductivity_vertical)
        self.Ss = self._float_or_none(specific_storage)
        self.Sy = self._float_or_none(specific_yield)
        self.well_top = self._float_or_none(well_screen_elevation_top)
        self.well_bot = self._float_or_none(well_screen_elevation_bottom)

        if (
            water_table_elevation is not None
            and saturated_thickness is not None
        ):
            raise RuntimeError(
                "RadialUnconfinedDrawdown() must specify only "
                + "water_table_elevation or saturated_thickness, but not "
                + "both at the same time."
            )

        if water_table_elevation is not None:
            self.saturated_thickness = (
                float(water_table_elevation) - self.bottom
            )
        elif saturated_thickness is not None:
            self.saturated_thickness = float(saturated_thickness)
        else:
            self.saturated_thickness = None

    def _prop_check(self):
        error = []
        if self.Kr is None:
            error.append("hydraulic_conductivity_radial")
        if self.Kz is None:
            error.append("hydraulic_conductivity_vertical")
        if self.Ss is None:
            error.append("specific_storage")
        if self.Sy is None:
            error.append("specific_yield")
        if self.well_top is None:
            error.append("well_screen_elevation_top")
        if self.well_bot is None:
            error.append("well_screen_elevation_bottom")
        if error:
            raise RuntimeError(
                "RadialUnconfinedDrawdown: Attempted to solve radial "
                + "groundwater model\nwith the following input not specified\n"
                + "\n".join(error)
            )
        if self.well_top <= self.well_bot:
            raise RuntimeError(
                "RadialUnconfinedDrawdown: "
                + "well_screen_elevation_top <= well_screen_elevation_bottom\n"
                + f"That is: {well_screen_elevation_top} <= "
                + f"{well_screen_elevation_bottom}"
            )

    def drawdown(
        self,
        pump,
        time,
        radius,
        observation_elevation,
        observation_elevation_bot=None,
        sumrtol=1.0e-6,
        u_n_rtol=1.0e-5,
        epsabs=1.49e-8,
        bessel_loop_limit=5,
        quad_limit=128,
        show_progress=False,
        ty_time=False,
        ts_time=False,
        as_head=False,
    ):
        """
        Solves the radial model's drawdown for a given pumping rate and
        time at a given observation point
        (radius, observation_elevation) or observation well screen interval
        (radius, observation_elevation:observation_elevation_bot).
        This solves drawdown by integrating equation 17 from
        Neuman, S. P. (1974). Effect of partial penetration on flow in
        unconfined aquifers considering delayed gravity response.
        Water resources research, 10(2), 303-312

        Parameters
        ----------
        pump : float
            Pumping rate of well at center of radial model ($L^3/T$)
            Positive values are the water extraction rate.
            Negative or zero values indicate no pumping and result returns
            the dimensionless drawdown instead of regular drawdown.
        time : float or Sequence[float]
            Time that observation is made
        radius : float
            Radius of the observation location (distance from well, $L$)
        observation_elevation : float
            Either the location of the observation point, or the top elevation
            of the observation well screen ($L$)
        observation_elevation_bot : float
            If specified, then represents the bottom elevation of the
            observation well screen. If not specified (or set to None), then
            observation location is treated as a single point, located at
            radius and observation_elevation ($L$)
        sumrtol : float
            Solution involves integration of $y$ variable from 0 to ∞ from
            Equation 17 in:
            Neuman, S. P. (1974). Effect of partial penetration on flow in
            unconfined aquifers considering delayed gravity response.
            Water resources research, 10(2), 303-312.

            The integration is broken into subsections that are spaced around
            bessel function roots. The integration is complete when a
            three sequential subsection solutions are less than
            sumrtol times the largest subsection.
            That is, the last included subsection contributes a
            relatively small value compared to the largest of the sum.
        u_n_rtol : float
            Terminates the solution of the infinite series:
            $\\sum_{n=1}^{\\infty} u_n(y)$
            when
            $u_n(y) < u_n(0) * u_n_rtol$
        epsabs : float or int
            scipy.integrate.quad absolute error tolerance.
            Passed directly to that function's `epsabs` kwarg.
        bessel_loop_limit : int
            the integral is solved along each bessel function root.
            The first 1024 roots are precalculated and automatically increased
            if more are required. The upper limit for calculated roots is
            1024 * 2 ^ bessel_loop_limit
            If this limit is reached, then a warning is raised.
        quad_limit : int
            scipy.integrate.quad upper bound on the number of
            subintervals used in the adaptive algorithm.
            Passed directly to that function's `limit` kwarg.
        show_progress : bool
           if True, then progress is printed to the command prompt in the form:
        ty_time : bool
           if True, then `time` kwarg is dimensionless time with
           respect to Specific Yield
        ts_time : bool
           if True, then `time` kwarg is dimensionless time with
           respect to Specific Storage.
        as_head : bool
            If true, then drawdown result is converted to
            head using the model bottom and initial saturated thickness.
            If pump > 0, then as_head is ignored.

        Returns
        -------
        result : float or list[float]
            If time is float, then result is float.
            If time is Sequence[float], then result is list[float].

            If pump > 0, then result is the drawdown that occurs
            from pump at time and radius at observation point
            observation_elevation or from the observation well
            screen interval observation_elevation to
            observation_elevation_top ($L$).


            If pump <= 0, then result is converted to
            dimensionless drawdown ($-$)
        """
        if not hasattr(time, "strip") and hasattr(time, "__iter__"):
            return self.drawdown_times(
                pump,
                time,
                radius,
                observation_elevation,
                observation_elevation_bot,
                sumrtol,
                u_n_rtol,
                epsabs,
                bessel_loop_limit,
                quad_limit,
                show_progress,
                ty_time,
                ts_time,
                as_head,
            )

        return self.drawdown_times(
            pump,
            [time],
            radius,
            observation_elevation,
            observation_elevation_bot,
            sumrtol,
            u_n_rtol,
            epsabs,
            bessel_loop_limit,
            quad_limit,
            show_progress,
            ty_time,
            ts_time,
            as_head,
        )[0]

    def drawdown_times(
        self,
        pump,
        times,
        radius,
        observation_elevation,
        observation_elevation_bot=None,
        sumrtol=1.0e-6,
        u_n_rtol=1.0e-5,
        epsabs=1.49e-8,
        bessel_loop_limit=5,
        quad_limit=128,
        show_progress=False,
        ty_time=False,
        ts_time=False,
        as_head=False,
    ):
        # Same as self.drawdown, but times is a list[float] of
        # observation times and returns a list[float] drawdowns.

        if bessel_loop_limit < 1:
            bessel_loop_limit = 1

        bessel_roots0 = 1024
        bessel_roots = bessel_roots0
        bessel_root_limit_reached = []

        self._prop_check()
        if ty_time and ts_time:
            raise RuntimeError(
                "RadialUnconfinedDrawdown.drawdown_times "
                + "cannot set both ty_time and ts_time to True."
            )

        r = radius
        b = self.saturated_thickness

        sigma = self.Ss * b / self.Sy
        beta = (r / b) * (r / b) * (self.Kz / self.Kr)
        sqrt_beta = sqrt(beta)

        if np.isnan(pump) or pump <= 0.0:
            # Return dimensionless drawdown
            coef = 1.0
        else:
            coef = pump / (4.0 * pi * b * self.Kr)

        # dimensionless well screen top
        dd = (self.saturated_thickness + self.bottom - self.well_top) / b
        # dimensionless well screen bottom
        ld = (self.saturated_thickness + self.bottom - self.well_bot) / b

        # Solution must be in dimensionless time with respect to Ss;
        # ts = kr*b*t/(Ss*b*r^2)
        if ty_time:
            ts_list = self.ty2ts(times)
        elif ts_time:
            ts_list = times
        else:
            ts_list = self.time2ts(times, r)

        # distance above bottom to observation point or obs screen bottom
        zt = observation_elevation - self.bottom
        if observation_elevation_bot is None:
            # Single Point Observation
            zd = zt / b  # dimensionless elevation of observation point
            neuman1974_integral = self.neuman1974_integral1
            obs_arg = (zd,)
        else:
            # distance above bottom to observation screen top
            zb = observation_elevation_bot - self.bottom
            # dimensionless elevation of observation screen interval
            ztd, zbd = zt / b, zb / b
            # dz = 1 / (zt - zb)  -> implied in the
            #                        modified u0 and uN functions
            neuman1974_integral = self.neuman1974_integral2
            obs_arg = (zbd, ztd)

        s = []  # drawdown, one to one match with times
        nstp = len(ts_list)
        for stp, ts in enumerate(ts_list):
            if show_progress:
                print(
                    f"Solving {stp+1:4d} of {nstp}; "
                    + f"time = {self.ts2time(ts, r)}",
                    end="",
                )

            args = (sigma, beta, sqrt_beta, ld, dd, ts, *obs_arg, u_n_rtol)
            sol = 0.0
            y0, y1 = 0.0, 0.0
            mxdelt = 0.0

            j0_roots = jn_zeros(0, bessel_roots) / sqrt_beta
            jr0 = 0
            jr1 = j0_roots.size

            converged = 0
            bessel_loop_count = 0
            while converged < 3 and bessel_loop_count <= bessel_loop_limit:
                if bessel_loop_count > 0:
                    bessel_roots *= 2
                    j0_roots = jn_zeros(0, bessel_roots) / sqrt_beta
                    jr0, jr1 = jr1, j0_roots.size

                j0_roots_iter = np.nditer(j0_roots[jr0:jr1])
                bessel_loop_count += 1
                # Iterate over two roots to get full cycle
                for j0_root in j0_roots_iter:
                    # First root
                    y0, y1 = y1, j0_root
                    delt1 = quad(
                        neuman1974_integral,
                        y0,
                        y1,
                        args,
                        epsabs=epsabs,
                        limit=quad_limit,
                    )[0]
                    #
                    # Second root
                    y0, y1 = y1, next(j0_roots_iter)
                    delt2 = quad(
                        neuman1974_integral,
                        y0,
                        y1,
                        args,
                        epsabs=epsabs,
                        limit=quad_limit,
                    )[0]

                    if np.isnan(delt1) or np.isnan(delt2):
                        break

                    sol += delt1 + delt2

                    adelt = abs(delt1 + delt2)
                    if adelt > mxdelt:
                        mxdelt = adelt
                    elif adelt < mxdelt * sumrtol:
                        converged += 1  # increment the convergence counter
                        # Converged if three sequential solutions (adelt)
                        # are less than mxdelt*sumrtol
                        if converged >= 3:
                            break
                    else:
                        converged = 0  # reset convergence counter
            if sol < 0.0:
                s.append(0.0)
            else:
                s.append(coef * sol)

            if converged < 3:
                bessel_root_limit_reached.append(stp)

            if show_progress:
                if converged < 3:
                    print(f"\ts = {s[-1]}\tbessel_loop_limit reached")
                else:
                    print(f"\ts = {s[-1]}")

        if pump > 0.0 and as_head:
            initial_head = self.bottom + self.saturated_thickness
            return [initial_head - drawdown for drawdown in s]

        if len(bessel_root_limit_reached) > 0:
            import warnings

            root = j0_roots[-1]
            bad_times = "\n".join(
                [str(times[it]) for it in bessel_root_limit_reached]
            )
            warnings.warn(
                f"\n\nRadialUnconfinedDrawdown.drawdown_times failed to "
                + f"meet convergence sumrtol = {sumrtol}"
                + "\nwithin the precalculated Bessel root solutions "
                + "(convergence is evaluated at every second Bessel root).\n\n"
                + "The number of Bessel roots are automatically increased "
                + "up to:\n"
                + f"   {bessel_roots0} * 2^bessel_loop_limit\nwhere:\n"
                + "   bessel_loop_limit = {bessel_loop_limit}\n"
                + f"resulting in {1024*2**bessel_loop_limit} roots evaluated, "
                + "with the last root being {root}\n"
                + f"(That is, the Neuman integral was solved form 0 to {root})"
                + "\n\n"
                + "You can either ignore this warning\n"
                + "or to remove it attempt to increase bessel_loop_limit\n"
                + "or increase sumrtol (reducing accuracy).\n\nThe following "
                + "times are what triggered this warning:\n"
                + bad_times
                + "\n"
            )
        return s

    @staticmethod
    def neuman1974_integral1(y, σ, β, sqrt_β, ld, dd, ts, zd, uN_tol=1.0e-6):
        """
        Solves equation 17 from
        Neuman, S. P. (1974). Effect of partial penetration on flow in
        unconfined aquifers considering delayed gravity response.
        Water resources research, 10(2), 303-312.
        """
        if y == 0.0 or ts == 0.0:
            return 0.0

        u0 = RadialUnconfinedDrawdown.u_0(σ, β, zd, ld, dd, ts, y)

        if np.isnan(u0):
            u0 = 0.0

        uN_func = RadialUnconfinedDrawdown.u_n
        mxdelt = 0.0
        uN = 0.0
        for n in range(1, 25001):
            delt = uN_func(σ, β, zd, ld, dd, ts, y, n)
            if np.isnan(delt):
                break
            uN += delt
            adelt = abs(delt)
            if adelt > mxdelt:
                mxdelt = adelt
            elif adelt < mxdelt * uN_tol:
                break

        return 4.0 * y * j0(y * sqrt_β) * (u0 + uN)

    @staticmethod
    def gamma0(g, y, s):
        """
        Gamma0 root function from equation 18 in:
        Neuman, S. P. (1974). Effect of partial penetration on flow in
        unconfined aquifers considering delayed gravity response.
        Water resources research, 10(2), 303-312.
            => Solution must be constrained by g^2 < y^2

        To honor the constraint solution returns the absolute value
         of the solution.
        """
        if g >= _hyperbolic_equivalence:
            # sinh ≈ cosh for large g
            return s * g - (y * y - g * g)

        return s * g * sinh(g) - (y * y - g * g) * cosh(g)

    @staticmethod
    def gammaN(g, y, s):
        """
        GammaN root function from equation 19 in:
        Neuman, S. P. (1974). Effect of partial penetration on flow in
        unconfined aquifers considering delayed gravity response.
        Water resources research, 10(2), 303-312.
            => Solution must be constrained by (2n-1)(π/2)< g < nπ
        """
        return s * g * sin(g) + (y * y + g * g) * cos(g)

    @staticmethod
    def u_0(σ, β, z, l, d, ts, y):
        gamma0 = RadialUnconfinedDrawdown.gamma0

        a, b = 0.9 * y, y
        try:
            a, b = RadialUnconfinedDrawdown._get_bracket(gamma0, a, b, (y, σ))
        except RuntimeError:
            a, b = RadialUnconfinedDrawdown._get_bracket(
                gamma0, 0.0, b, (y, σ), 1000
            )

        g = brentq(gamma0, a, b, args=(y, σ), maxiter=500, xtol=1.0e-16)

        # Check for cosh/sinh overflow
        if g > _hyperbolic_max_value:
            return 0.0

        y2 = y * y
        g2 = g * g
        num1 = 1 - exp(-ts * β * (y2 - g2))
        num2 = cosh(g * z)
        num3 = sinh(g * (1 - d)) - sinh(g * (1 - l))
        den1 = y2 + (1 + σ) * g2 - ((y2 - g2) ** 2) / σ
        den2 = cosh(g)
        den3 = (l - d) * sinh(g)
        # num1*num2*num3 / (den1*den2*den3)
        return (num1 / den1) * (num2 / den2) * (num3 / den3)

    @staticmethod
    def u_n(σ, β, z, l, d, ts, y, n):
        gammaN = RadialUnconfinedDrawdown.gammaN

        a, b = (2 * n - 1) * (pi / 2.0), n * pi
        try:
            a, b = RadialUnconfinedDrawdown._get_bracket(gammaN, a, b, (y, σ))
        except RuntimeError:
            a, b = RadialUnconfinedDrawdown._get_bracket(
                gammaN, a, b, (y, σ), 1000
            )

        g = brentq(gammaN, a, b, args=(y, σ), maxiter=500, xtol=1.0e-16)

        y2 = y * y
        g2 = g * g
        num1 = 1 - exp(-ts * β * (y2 + g2))
        num2 = cos(g * z)
        num3 = sin(g * (1 - d)) - sin(g * (1 - l))
        den1 = y2 - (1 + σ) * g2 - ((y2 + g2) ** 2) / σ
        den2 = cos(g)
        den3 = (l - d) * sin(g)
        return num1 * num2 * num3 / (den1 * den2 * den3)

    @staticmethod
    def neuman1974_integral2(
        y, σ, β, sqrt_β, ld, dd, ts, z1, z2, uN_tol=1.0e-10
    ):
        """
        Solves equation 20 from
        Neuman, S. P. (1974). Effect of partial penetration on flow in
        unconfined aquifers considering delayed gravity response.
        Water resources research, 10(2), 303-312.
        """
        if y == 0.0 or ts == 0.0:
            return 0.0

        u0 = RadialUnconfinedDrawdown.u_0_z1z2(σ, β, z1, z2, ld, dd, ts, y)

        uN_func = RadialUnconfinedDrawdown.u_n_z1z2
        mxdelt = 0.0
        uN = 0.0
        for n in range(1, 10001):
            delt = uN_func(σ, β, z1, z2, ld, dd, ts, y, n)
            uN += delt
            adelt = abs(delt)
            if adelt > mxdelt:
                mxdelt = adelt
            elif adelt < mxdelt * uN_tol:
                break

        return 4.0 * y * j0(y * sqrt_β) * (u0 + uN)

    @staticmethod
    def u_0_z1z2(σ, β, z1, z2, l, d, ts, y):
        gamma0 = RadialUnconfinedDrawdown.gamma0

        a, b = 0.9 * y, y
        try:
            a, b = RadialUnconfinedDrawdown._get_bracket(gamma0, a, b, (y, σ))
        except RuntimeError:
            a, b = RadialUnconfinedDrawdown._get_bracket(
                gamma0, 0.0, b, (y, σ), 1000
            )

        g = brentq(gamma0, a, b, args=(y, σ), maxiter=500, xtol=1.0e-16)

        # Check for cosh/sinh overflow
        if g > _hyperbolic_max_value:
            return 0.0

        y2 = y * y
        g2 = g * g
        num1 = 1 - exp(-ts * β * (y2 - g2))
        num2 = sinh(g * z2) - sinh(g * z1)
        num3 = sinh(g * (1 - d)) - sinh(g * (1 - l))
        den1 = (y2 + (1 + σ) * g2 - ((y2 - g2) ** 2) / σ) * (z2 - z1) * g
        den2 = cosh(g)
        den3 = (l - d) * sinh(g)
        # num1*num2*num3 / (den1*den2*den3)
        return (num1 / den1) * (num2 / den2) * (num3 / den3)

    @staticmethod
    def u_n_z1z2(σ, β, z1, z2, l, d, ts, y, n):
        gammaN = RadialUnconfinedDrawdown.gammaN

        a, b = (2 * n - 1) * (pi / 2.0), n * pi
        try:
            a, b = RadialUnconfinedDrawdown._get_bracket(gammaN, a, b, (y, σ))
        except RuntimeError:
            a, b = RadialUnconfinedDrawdown._get_bracket(
                gammaN, a, b, (y, σ), 1000
            )

        g = brentq(gammaN, a, b, args=(y, σ), maxiter=500, xtol=1.0e-16)

        y2 = y * y
        g2 = g * g
        num1 = 1 - exp(-ts * β * (y2 + g2))
        num2 = sin(g * z2) - sin(g * z1)
        num3 = sin(g * (1 - d)) - sin(g * (1 - l))
        den1 = y2 - (1 + σ) * g2 - ((y2 + g2) ** 2) / σ
        den2 = cos(g) * (z2 - z1) * g
        den3 = (l - d) * sin(g)
        return num1 * num2 * num3 / (den1 * den2 * den3)

    def time2ty(self, time, radius):
        # dimensionless time with respect to Sy
        if hasattr(time, "__iter__"):
            # can iterate to get multiple times
            return [
                self.Kr
                * self.saturated_thickness
                * t
                / (self.Sy * radius * radius)
                for t in time
            ]
        return (
            self.Kr
            * self.saturated_thickness
            * time
            / (self.Sy * radius * radius)
        )

    def time2ts(self, time, radius):
        # dimensionless time with respect to Ss
        if hasattr(time, "__iter__"):
            # can iterate to get multiple times
            return [self.Kr * t / (self.Ss * radius * radius) for t in time]
        return self.Kr * time / (self.Ss * radius * radius)

    def ty2time(self, ty, radius):
        # dimensionless time with respect to Sy
        if hasattr(ty, "__iter__"):
            # can iterate to get multiple times
            return [
                t
                * self.Sy
                * radius
                * radius
                / (self.Kr * self.saturated_thickness)
                for t in ty
            ]
        return (
            ty
            * self.Sy
            * radius
            * radius
            / (self.Kr * self.saturated_thickness)
        )

    def ts2time(self, ts, radius):  # dimensionless time with respect to Ss
        if hasattr(ts, "__iter__"):  # can iterate to get multiple times
            return [t * self.Ss * radius * radius / self.Kr for t in ts]
        return ts * self.Ss * radius * radius / self.Kr

    def ty2ts(self, ty):
        if hasattr(ty, "__iter__"):
            # can iterate to get multiple times
            return [
                t * self.Sy / (self.Ss * self.saturated_thickness) for t in ty
            ]
        return ty * self.Sy / (self.Ss * self.saturated_thickness)

    def drawdown2unitless(self, s, pump):
        # dimensionless drawdown
        return 4 * pi * self.Kr * self.saturated_thickness * s / pump

    def unitless2drawdown(self, s, pump):
        # drawdown
        return pump * s / (4 * pi * self.Kr * self.saturated_thickness)

    @staticmethod
    def _float_or_none(val):
        if val is not None:
            return float(val)
        return None

    @staticmethod
    def _get_bracket(func, a, b, arg=(), internal_search_split=100):
        """
        Given initial range [a, b], search within the range for
        root finding brackets.
        That is, return [a, b] that results in f(a) * f(b) < 0.
        """
        if a > b:
            a, b = b, a

        f1 = func(a, *arg)
        f2 = func(b, *arg)

        if f1 * f2 <= 0.0:
            return a, b

        # same sign, search within for sign change
        delt = abs(b - a) / internal_search_split
        a -= delt
        for _ in range(internal_search_split):
            a += delt
            f1 = func(a, *arg)
            if f1 * f2 <= 0.0:
                return a, b

        raise RuntimeError(
            "get_bracket: failed to find bracket interval with opposite "
            + f"signs, that is: f(a)*f(b) < 0 for func: {func}"
        )


if __name__ == "__main__":
    # Example validation using Figure 2 from
    # Neuman, S. P. (1974). Effect of partial penetration on flow in
    #     unconfined aquifers considering delayed gravity response.
    #     Water resources research, 10(2), 303-312.
    # (no units are provided because result is dimensionless)

    b = 100  # saturated thickness
    bottom = 0.0
    well_top = b + bottom  # Top of well is top of aquifer
    well_bot = (bottom + b) - 0.2 * b  # Well bottom is 20% below aquifer top
    ss = 1.0e-05
    sy = 0.1
    kz = 10.0
    kr = kz
    kd = kz / kr
    r = 0.6 * b / (kd**0.5)
    z = 0.85 * b
    # beta = (r / b) * (r / b) * (kz / kr)  # = .36
    pump = 1000.0

    Neuman1974 = RadialUnconfinedDrawdown(
        bottom_elevation=bottom,
        hydraulic_conductivity_radial=kr,
        hydraulic_conductivity_vertical=kz,
        specific_storage=ss,
        specific_yield=sy,
        well_screen_elevation_top=well_top,
        well_screen_elevation_bottom=well_bot,
        saturated_thickness=b,
    )

    # Convert ts to regular time
    times = [Neuman1974.ts2time(ts, r) for ts in [1.0, 10.0, 100.0]]

    s1 = Neuman1974.drawdown(pump, times[0], r, z)  # ≈ 0.0092
    s2 = Neuman1974.drawdown(pump, times[1], r, z)  # ≈ 0.0156
    s3 = Neuman1974.drawdown(pump, times[2], r, z)  # ≈ 0.0733
    assert int(s1 * 10000) / 10000 == 0.0092
    assert int(s2 * 10000) / 10000 == 0.0156
    assert int(s3 * 10000) / 10000 == 0.0733

    sd1 = Neuman1974.drawdown_dimensionless(times[0], r, z)  # ≈ 0.116
    sd2 = Neuman1974.drawdown_dimensionless(times[1], r, z)  # ≈ 0.196
    sd3 = Neuman1974.drawdown_dimensionless(times[2], r, z)  # ≈ 0.921
    assert int(sd1 * 1000) / 1000 == 0.116
    assert int(sd2 * 1000) / 1000 == 0.196
    assert int(sd3 * 1000) / 1000 == 0.921

    ts_time = np.logspace(-2, 3, 20)
    times = [
        Neuman1974.ts2time(ts, r) for ts in ts_time
    ]  # Convert ts to regular time
    #
    # pump < 0 ⇒ dimensionless drawdown
    sd = Neuman1974.drawdown(
        -1.0, ts_time, r, z, ts_time=True, show_progress=True
    )
    s = [pump * dd / (4 * pi * kr * b) for dd in sd]

    # From digitizing Figure 2 from Neuman 1974
    ts_true = [
        0.089687632,
        0.091225896,
        0.093981606,
        0.099322056,
        0.105413236,
        0.110460061,
        0.116736935,
        0.123370444,
        0.134891669,
        0.149382,
        0.166841577,
        0.18087786,
        0.197769543,
        0.213498224,
        0.236432628,
        0.264066571,
        0.292433127,
        0.322472937,
        0.363238602,
        0.403972132,
        0.466797384,
        0.541691087,
        0.639382207,
        0.804387135,
        0.994911098,
        1.300487774,
        1.561350982,
        1.796516091,
        2.289153019,
        2.929308713,
        3.623130394,
        4.61665772,
        5.637774401,
        7.338096581,
        8.735440623,
        10.57723024,
        12.37907737,
        14.18313516,
        17.24669923,
        19.67631134,
        23.32377322,
        25.82925931,
        29.59349474,
        33.47660549,
        39.34627872,
        45.46538164,
        53.21043662,
        58.92640751,
        70.14741888,
        80.02936879,
        93.66241967,
        107.7695188,
        122.4298652,
        139.084521,
        158.0047806,
        182.5775622,
        214.590086,
        279.309292,
        362.0049489,
        469.1848385,
        615.9033614,
        808.5008911,
        1000.0,
    ]
    sd_true = [
        0.010442884,
        0.011437831,
        0.012691503,
        0.014834303,
        0.017114926,
        0.019239359,
        0.021627491,
        0.024630197,
        0.028664212,
        0.0336493,
        0.040018293,
        0.044404599,
        0.049271676,
        0.053966039,
        0.060141169,
        0.065302682,
        0.071524401,
        0.075998143,
        0.082163706,
        0.087302955,
        0.094795544,
        0.1007249,
        0.107489912,
        0.114709285,
        0.121884175,
        0.130635227,
        0.134076495,
        0.138206013,
        0.142462775,
        0.150067308,
        0.156035929,
        0.164364984,
        0.173890601,
        0.187997964,
        0.202370903,
        0.214099117,
        0.232473715,
        0.245946509,
        0.274089351,
        0.296325619,
        0.323154565,
        0.350888753,
        0.376081747,
        0.419115531,
        0.475240418,
        0.529620599,
        0.585129692,
        0.654916226,
        0.749082106,
        0.849394041,
        0.958973727,
        1.073347435,
        1.206579326,
        1.32153628,
        1.460045395,
        1.620077321,
        1.805457382,
        2.156510354,
        2.488051985,
        2.833486709,
        3.226880699,
        3.580571026,
        3.804531016,
    ]
    fig = plt.figure(figsize=(6, 4))
    ax = fig.add_subplot()
    ax.set_title("Neuman 1974 Figure 2 from Python")
    ax.set_xlabel("$t_s$ (-)")
    ax.set_ylabel("$s_d$ (-)")
    # plt.loglog(ts_true, sd_true)
    plt.loglog(
        ts_true, sd_true, "-k", markerfacecolor="none", label="Neuman Solution"
    )
    plt.loglog(
        ts_time, sd, "ob", markerfacecolor="none", label="Python Solution"
    )
    ax.set_xlim(0.01, 1000.0)
    ax.set_ylim(0.01, 10.0)
    ax.legend()  # loc=2
    ax.grid(visible=True, which="major", axis="both")
    ax.set_zorder(2)
    ax.set_facecolor("none")

    plt.tight_layout()
    plt.show()
    pass
