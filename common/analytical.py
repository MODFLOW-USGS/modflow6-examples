import numpy as np
from scipy.special import erfc
from scipy.special import erf


def diffusion(x, t, v, R, D):
    """
    Calculate the analytical solution for one-dimension advection and
    dispersion using the solution of Lapidus and Amundson (1952) and
    Ogata and Banks (1961)

    Parameters
    ----------
    x : float or ndarray
        x position
    t : float or ndarray
        time
    v : float or ndarray
        velocity
    R : float or ndarray
        retardation factor, a value of one indicates there is no adsorption
    D : float or ndarray
        diffusion coefficient in length squared per time

    Returns
    -------
    result : float or ndarray
        normalized concentration value

    """
    denom = 2.0 * np.sqrt(D * R * t)
    t1 = 0.5 * erfc((R * x - v * t) / denom)
    t2 = 0.5 * np.exp(v * x / D)
    t3 = erfc((R * x + v * t) / denom)
    return t1 + t2 * t3


class Wexler1d(object):
    """
    Analytical solution for 1D transport with inflow at a concentration of 1.
    at x=0 and a third-type bound at location l.
    Wexler Page 17 and Van Genuchten and Alves pages 66-67
    """

    def betaeqn(self, beta, d, v, l):
        return beta / np.tan(beta) - beta ** 2 * d / v / l + v * l / 4.0 / d

    def fprimebetaeqn(self, beta, d, v, l):
        """
        f1 = cotx - x/sinx2 - (2.0D0*C*x)

        """
        c = v * l / 4.0 / d
        return 1.0 / np.tan(beta) - beta / np.sin(beta) ** 2 - 2.0 * c * beta

    def fprime2betaeqn(self, beta, d, v, l):
        """
        f2 = -1.0D0/sinx2 - (sinx2-x*DSIN(x*2.0D0))/(sinx2*sinx2) - 2.0D0*C

        """
        c = v * l / 4.0 / d
        sinx2 = np.sin(beta) ** 2
        return (
            -1.0 / sinx2
            - (sinx2 - beta * np.sin(beta * 2.0)) / (sinx2 * sinx2)
            - 2.0 * c
        )

    def solvebetaeqn(self, beta, d, v, l, xtol=1.0e-12):
        from scipy.optimize import fsolve

        t = fsolve(
            self.betaeqn,
            beta,
            args=(d, v, l),
            fprime=self.fprime2betaeqn,
            xtol=xtol,
            full_output=True,
        )
        result = t[0][0]
        infod = t[1]
        isoln = t[2]
        msg = t[3]
        if abs(result - beta) > np.pi:
            raise Exception("Error in beta solution")
        err = self.betaeqn(result, d, v, l)
        fvec = infod["fvec"][0]
        if isoln != 1:
            print("Error in beta solve", err, result, d, v, l, msg)
        return result

    def root3(self, d, v, l, nval=1000):
        b = 0.5 * np.pi
        betalist = []
        for i in range(nval):
            b = self.solvebetaeqn(b, d, v, l)
            err = self.betaeqn(b, d, v, l)
            betalist.append(b)
            b += np.pi
        return betalist

    def analytical(self, x, t, v, l, d, tol=1.0e-20, nval=5000):
        sigma = 0.0
        betalist = self.root3(d, v, l, nval=nval)
        for i, bi in enumerate(betalist):

            denom = bi ** 2 + (v * l / 2.0 / d) ** 2 + v * l / d
            x1 = (
                bi
                * (
                    bi * np.cos(bi * x / l)
                    + v * l / 2.0 / d * np.sin(bi * x / l)
                )
                / denom
            )

            denom = bi ** 2 + (v * l / 2.0 / d) ** 2
            x2 = np.exp(-1 * bi ** 2 * d * t / l ** 2) / denom

            sigma += x1 * x2
            term1 = (
                2.0
                * v
                * l
                / d
                * np.exp(v * x / 2.0 / d - v ** 2 * t / 4.0 / d)
            )
            conc = 1.0 - term1 * sigma
            if i > 0:
                diff = abs(conc - concold)
                if np.all(diff < tol):
                    break
            concold = conc
        return conc

    def analytical2(self, x, t, v, l, d, e=0.0, tol=1.0e-20, nval=5000):
        """
        Calculate the analytical solution for one-dimension advection and
        dispersion using the solution of Lapidus and Amundson (1952) and
        Ogata and Banks (1961)

        Parameters
        ----------
        x : float or ndarray
            x position
        t : float or ndarray
            time
        v : float or ndarray
            velocity
        l : float
            length domain
        d : float
            dispersion coefficient
        e : float
            decay rate

        Returns
        -------
        result : float or ndarray
            normalized concentration value

        """
        u = v ** 2 + 4.0 * e * d
        u = np.sqrt(u)
        sigma = 0.0
        denom = (u + v) / 2.0 / v - (u - v) ** 2.0 / 2.0 / v / (
            u + v
        ) * np.exp(-u * l / d)
        term1 = np.exp((v - u) * x / 2.0 / d) + (u - v) / (u + v) * np.exp(
            (v + u) * x / 2.0 / d - u * l / d
        )
        term1 = term1 / denom
        term2 = (
            2.0
            * v
            * l
            / d
            * np.exp(v * x / 2.0 / d - v ** 2 * t / 4.0 / d - e * t)
        )
        betalist = self.root3(d, v, l, nval=nval)
        for i, bi in enumerate(betalist):

            denom = bi ** 2 + (v * l / 2.0 / d) ** 2 + v * l / d
            x1 = (
                bi
                * (
                    bi * np.cos(bi * x / l)
                    + v * l / 2.0 / d * np.sin(bi * x / l)
                )
                / denom
            )

            denom = bi ** 2 + (v * l / 2.0 / d) ** 2 + e * l ** 2 / d
            x2 = np.exp(-1 * bi ** 2 * d * t / l ** 2) / denom

            sigma += x1 * x2

            conc = term1 - term2 * sigma
            if i > 0:
                diff = abs(conc - concold)
                if np.all(diff < tol):
                    break
            concold = conc
        return conc


class Wexler3d:
    """
    Analytical solution for 3D transport with inflow at a well with a
    specified concentration.
    Wexler Page 47
    """

    def calcgamma(self, x, y, z, xc, yc, zc, dx, dy, dz):
        gam = np.sqrt(
            (x - xc) ** 2 + dx / dy * (y - yc) ** 2 + dx / dz * (z - zc) ** 2
        )
        return gam

    def calcbeta(self, v, dx, gam, lam):
        beta = np.sqrt(v ** 2 + 4.0 * dx * gam * lam)
        return beta

    def analytical(
        self, x, y, z, t, v, xc, yc, zc, dx, dy, dz, n, q, lam=0.0, c0=1.0
    ):
        gam = self.calcgamma(x, y, z, xc, yc, zc, dx, dy, dz)
        beta = self.calcbeta(v, dx, gam, lam)
        term1 = (
            c0
            * q
            * np.exp(v * (x - xc) / 2.0 / dx)
            / 8.0
            / n
            / np.pi
            / gam
            / np.sqrt(dy * dz)
        )
        term2 = np.exp(gam * beta / 2.0 / dx) * erfc(
            (gam + beta * t) / 2.0 / np.sqrt(dx * t)
        )
        term3 = np.exp(-gam * beta / 2.0 / dx) * erfc(
            (gam - beta * t) / 2.0 / np.sqrt(dx * t)
        )
        return term1 * (term2 + term3)

    def multiwell(
        self, x, y, z, t, v, xc, yc, zc, dx, dy, dz, n, ql, lam=0.0, c0=1.0
    ):
        shape = self.analytical(
            x, y, z, t, v, xc[0], yc[0], zc[0], dx, dy, dz, n, ql[0], lam
        ).shape
        result = np.zeros(shape)
        for xx, yy, zz, q in zip(xc, yc, zc, ql):
            result += self.analytical(
                x, y, z, t, v, xx, yy, zz, dx, dy, dz, n, q, lam, c0
            )
        return result


class BakkerRotatingInterface:
    """
    Analytical solution for rotating interfaces (Bakker et al. 2004)

    """

    @staticmethod
    def get_s(k, rhoa, rhob, alpha):
        return k * (rhob - rhoa) / rhoa * np.sin(alpha)

    @staticmethod
    def get_F(z, zeta1, omega1, s):
        l = (zeta1.real - omega1.real) ** 2 + (zeta1.imag - omega1.imag) ** 2
        l = np.sqrt(l)
        try:
            v = (
                s
                * l
                * complex(0, 1)
                / 2
                / np.pi
                / (zeta1 - omega1)
                * np.log((z - zeta1) / (z - omega1))
            )
        except:
            v = 0.0
        return v

    @staticmethod
    def get_Fgrid(xg, yg, zeta1, omega1, s):
        qxg = []
        qyg = []
        for x, y in zip(xg.flatten(), yg.flatten()):
            z = complex(x, y)
            W = BakkerRotatingInterface.get_F(z, zeta1, omega1, s)
            qx = W.real
            qy = -W.imag
            qxg.append(qx)
            qyg.append(qy)
        qxg = np.array(qxg)
        qyg = np.array(qyg)
        qxg = qxg.reshape(xg.shape)
        qyg = qyg.reshape(yg.shape)
        return qxg, qyg

    @staticmethod
    def get_zetan(n, x0, a, b):
        return complex(x0 + (-1) ** n * a, (2 * n - 1) * b)

    @staticmethod
    def get_omegan(n, x0, a, b):
        return complex(x0 + (-1) ** (1 + n) * a, -(2 * n - 1) * b)

    @staticmethod
    def get_w(xg, yg, k, rhoa, rhob, a, b, x0):
        zeta1 = BakkerRotatingInterface.get_zetan(1, x0, a, b)
        omega1 = BakkerRotatingInterface.get_omegan(1, x0, a, b)
        alpha = np.arctan2(b, a)
        s = BakkerRotatingInterface.get_s(k, rhoa, rhob, alpha)
        qxg, qyg = BakkerRotatingInterface.get_Fgrid(xg, yg, zeta1, omega1, s)
        for n in range(1, 5):
            zetan = BakkerRotatingInterface.get_zetan(n, x0, a, b)
            zetanp1 = BakkerRotatingInterface.get_zetan(n + 1, x0, a, b)
            qx1, qy1 = BakkerRotatingInterface.get_Fgrid(
                xg, yg, zetan, zetanp1, (-1) ** n * s
            )
            omegan = BakkerRotatingInterface.get_omegan(n, x0, a, b)
            omeganp1 = BakkerRotatingInterface.get_omegan(n + 1, x0, a, b)
            qx2, qy2 = BakkerRotatingInterface.get_Fgrid(
                xg, yg, omegan, omeganp1, (-1) ** n * s
            )
            qxg += qx1 + qx2
            qyg += qy1 + qy2
        return qxg, qyg


def diffusion(x, t, v, R, D):
    """
    Calculate the analytical solution for one-dimension advection and
    dispersion using the solution of Lapidus and Amundson (1952) and
    Ogata and Banks (1961)

    Parameters
    ----------
    x : float or ndarray
        x position
    t : float or ndarray
        time
    v : float or ndarray
        velocity
    R : float or ndarray
        retardation factor, a value of one indicates there is no adsorption
    D : float or ndarray
        diffusion coefficient in length squared per time

    Returns
    -------
    result : float or ndarray
        normalized concentration value

    """
    denom = 2.0 * np.sqrt(D * R * t)
    t1 = 0.5 * erfc((R * x - v * t) / denom)
    t2 = 0.5 * np.exp(v * x / D)
    t3 = erfc((R * x + v * t) / denom)
    return t1 + t2 * t3


def hechtMendez_SS_3d(
    x_pos, To, Y3d, Z3d, ath, atv, Fplanar, va, n, rhow, cw, thermdiff
):
    """
    Calculate the analytical solution for changes in temperature three-
    dimensional changes in temperature using transient solution provided in
    the appendix of Hecht-Mendez et al. (2010) as equation A5.  Note that for
    SS conditions, the erfc term reduces to 1 as t -> infinity and the To/2
    term becomes T.

    Parameters
    ----------
    x_pos : float or ndarray
        x position
    To : float or ndarray
         initial temperature of the ground, degrees K
    Y3d : float or ndarray
          dimension of source in y direction for 3D test problem
    Z3d : float or ndarray
          dimension of source in z direction for 3D test problem
    ath : float or ndarray
          transverse horizontal dispersivity
    atv : float or ndarray
          transverse vertical dispersivity
    Fplanar : float or ndarray
              energy extraction (point source)
    va : float or ndarray
         seepage velocity
    n : float or ndarray
        porosity
    rhow : float or ndarray
           desity of water
    cw : float or ndarray
         specific heat capacity of water
    thermdiff : float or ndarray
                molecular diffusion coefficient, or in this case thermal
                diffusivity
    """

    # calculate transverse horizontal heat dispersion
    Dy = ath * (va ** 2 / abs(va)) + thermdiff
    t2 = erf(Y3d / (4 * np.sqrt(Dy * (x_pos / va))))

    Dz = atv * (va ** 2 / abs(va)) + thermdiff
    t3 = erf(Z3d / (4 * np.sqrt(Dz * (x_pos / va))))

    # initial temperature at the source
    To_planar = Fplanar / (abs(va) * n * rhow * cw)

    sln = To + (To_planar * t2 * t3)
    return sln


def hechtMendezSS(x_pos, y, a, F0, va, n, rhow, cw, thermdiff):
    """
    Calculate the analytical solution for changes in temperature three-
    dimensional changes in temperature for a steady state solution provided in
    the appendix of Hecht-Mendez et al. (2010) as equation A4

    Parameters
    ----------
    x : float or ndarray
        x position
    y : float or ndarray
        y position
    a : float or ndarray
        longitudinal dispersivity
    F0 : float or ndarray
         energy extraction (point source)
    va : float or ndarray
         seepage velocity
    n : float or ndarray
        porosity
    rhow : float or ndarray
           desity of water
    cw : float or ndarray
         specific heat capacity of water
    thermdiff : float or ndarray
                molecular diffusion coefficient, or in this case thermal
                diffusivity
    """

    # calculate transverse horizontal heat dispersion
    Dth = a * (va ** 2 / abs(va)) + thermdiff

    t1 = F0 / (
        va * n * rhow * cw * ((4 * np.pi * Dth * (x_pos / va)) ** (0.5))
    )
    t2 = np.exp((-1 * va * y ** 2) / (4 * Dth * x_pos))
    sln = t1 * t2
    return sln


def hechtMendez3d(
    x_pos, t, Y, Z, al, ath, atv, thermdiff, va, n, R, Fplanar, cw, rhow
):
    """
    Calculate the analytical solution for three-dimensional changes in
    temperature based on the solution provided in the appendix of Hecht-Mendez
    et al. (2010) as equation A5

    Parameters
    ----------
    x : float or ndarray
        x position
    t : float or ndarray
        time
    Y : float or ndarray
        dimension of the source in the y direction
    Z : float or ndarray
        dimension of the source in the z direction
    al : float or ndarray
         longitudinal dispersivity
    ath : float or ndarray
          transverse horizontal dispersivity
    atv : float or ndarray
          transverse vertical dispersivity
    thermdiff : float or ndarray
                molecular diffusion coefficient, or in this case thermal
                diffusivity
    va : float or ndarray
         seepage velocity
    n : float or ndarray
        porosity
    R : float or ndarray
        retardation coefficient
    Fplanar : float or ndarray
              energy extraction (point source)
    cw : float or ndarray
         specific heat capacity of water
    rhow : float or ndarray
           desity of water

    """
    To_planar = Fplanar / (va * n * rhow * cw)

    Dl = al * (va ** 2 / abs(va)) + thermdiff
    numer = R * x_pos - va * t
    denom = 2 * np.sqrt(Dl * R * t)

    t1 = (To_planar / 2) * erfc(numer / denom)

    Dth = ath * (va ** 2 / abs(va)) + thermdiff
    t2 = erf(Y / (4 * np.sqrt(Dth * (x_pos / va))))

    Dtv = atv * (va ** 2 / abs(va)) + thermdiff
    t3 = erf(Z / (4 * np.sqrt(Dtv * (x_pos / va))))

    sln = t1 * t2 * t3
    return sln

