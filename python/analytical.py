import numpy as np
from scipy.special import erfc

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
    denom = 2. * np.sqrt(D * R * t)
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
        return beta / np.tan(beta) - beta ** 2 * d / v / l + v * l / 4. / d

    def fprimebetaeqn(self, beta, d, v, l):
        """
        f1 = cotx - x/sinx2 - (2.0D0*C*x)

        """
        c = v * l / 4. / d
        return 1. / np.tan(beta) - beta / np.sin(beta) ** 2 - 2. * c * beta

    def fprime2betaeqn(self, beta, d, v, l):
        """
        f2 = -1.0D0/sinx2 - (sinx2-x*DSIN(x*2.0D0))/(sinx2*sinx2) - 2.0D0*C

        """
        c = v * l / 4. / d
        sinx2 = np.sin(beta) ** 2
        return -1. / sinx2 - (sinx2 - beta * np.sin(beta * 2.)) / (sinx2 * sinx2) - 2. * c

    def solvebetaeqn(self, beta, d, v, l, xtol = 1.e-12):
        from scipy.optimize import fsolve
        t = fsolve(self.betaeqn, beta, args=(d, v, l),
                   fprime=self.fprime2betaeqn,
                   xtol=xtol, full_output=True)
        result = t[0][0]
        infod = t[1]
        isoln = t[2]
        msg = t[3]
        if abs(result - beta) > np.pi:
            raise Exception('Error in beta solution')
        err = self.betaeqn(result, d, v, l)
        fvec = infod['fvec'][0]
        if isoln != 1:
            print('Error in beta solve', err, result, d, v, l, msg)
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

    def analytical(self, x, t, v, l, d, tol=1.e-20, nval=5000):
        sigma = 0.
        betalist = self.root3(d, v, l, nval=nval)
        for i, bi in enumerate(betalist):
            
            denom = (bi ** 2 + (v * l / 2. / d) ** 2 + v * l / d)
            x1 = bi * (bi * np.cos(bi * x / l) + v * l / 2. / d * np.sin(bi * x / l)) / denom
            
            denom = (bi ** 2 + (v * l / 2. / d) ** 2)
            x2 = np.exp(-1 * bi ** 2 * d * t / l ** 2) / denom
            
            sigma += x1 * x2
            term1 = 2. * v * l / d * np.exp(v * x / 2. / d - v ** 2 * t / 4. / d)
            conc = 1. - term1 * sigma
            if i > 0:
                diff = abs(conc - concold)
                if np.all(diff < tol):
                    break
            concold = conc
        return conc
        
    def analytical2(self, x, t, v, l, d, e=0., tol=1.e-20, nval=5000):
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
        e : float
            decay rate
    
        Returns
        -------
        result : float or ndarray
            normalized concentration value
    
        """
        u = v ** 2 + 4. * e * d
        u = np.sqrt(u)
        sigma = 0.
        denom = (u + v) / 2. / v - (u - v) ** 2. / 2. / v / (u + v) * np.exp(-u * l / d)
        term1 = np.exp( (v - u) * x / 2. / d) + (u - v) / (u + v) * np.exp((v + u) * x / 2. / d - u * l / d)
        term1 = term1 / denom
        term2 = 2. * v * l / d * np.exp(v * x / 2. / d - v ** 2 * t / 4. / d - e * t)
        betalist = self.root3(d, v, l, nval=nval)
        for i, bi in enumerate(betalist):

            denom = (bi ** 2 + (v * l / 2. / d) ** 2 + v * l / d)
            x1 = bi * (bi * np.cos(bi * x / l) + v * l / 2. / d * np.sin(bi * x / l)) / denom
            
            denom = bi ** 2 + (v * l / 2. / d) ** 2 + e * l ** 2 / d
            x2 = np.exp(-1 * bi ** 2 * d * t / l ** 2) / denom
            
            sigma += x1 * x2
            
            conc = term1 - term2 * sigma
            if i > 0:
                diff = abs(conc - concold)
                if np.all(diff < tol):
                    break
            concold = conc
        return conc
        

class Wexler3d():
    """
    Analytical solution for 3D transport with inflow at a well with a
    specified concentration.
    Wexler Page 47
    """

    def calcgamma(self, x, y, z, xc, yc, zc, dx, dy, dz):
        gam = np.sqrt(
            (x - xc) ** 2 + dx / dy * (y - yc) ** 2 + dx / dz * (z - zc) ** 2)
        return gam

    def calcbeta(self, v, dx, gam, lam):
        beta = np.sqrt(v ** 2 + 4. * dx * gam * lam)
        return beta

    def analytical(self, x, y, z, t, v, xc, yc, zc, dx, dy, dz, n, q, lam=0., c0=1.):
        gam = self.calcgamma(x, y, z, xc, yc, zc, dx, dy, dz)
        beta = self.calcbeta(v, dx, gam, lam)
        term1 = c0 * q * np.exp(
            v * (x - xc) / 2. / dx) / 8. / n / np.pi / gam / np.sqrt(dy * dz)
        term2 = np.exp(gam * beta / 2. / dx) * erfc(
            (gam + beta * t) / 2. / np.sqrt(dx * t))
        term3 = np.exp(-gam * beta / 2. / dx) * erfc(
            (gam - beta * t) / 2. / np.sqrt(dx * t))
        return term1 * (term2 + term3)

    def multiwell(self, x, y, z, t, v, xc, yc, zc, dx, dy, dz, n, ql, lam=0., c0=1.):
        shape = self.analytical(x, y, z, t, v, xc[0], yc[0], zc[0], dx, dy, dz, n,
                           ql[0], lam).shape
        result = np.zeros(shape)
        for xx, yy, zz, q in zip(xc, yc, zc, ql):
            result += self.analytical(x, y, z, t, v, xx, yy, zz, dx, dy, dz, n, q,
                                 lam, c0)
        return result
