#!/usr/bin/env python
"""
Tests for measuring things

Run with:
   python KronPhotometry.py
or
   python
   >>> import KronPhotometry; KronPhotometry.run()
"""

import math, os, sys, unittest
import lsst.utils.tests as tests
import lsst.pex.exceptions as pexExceptions
import lsst.pex.logging as pexLogging
import lsst.pex.policy as pexPolicy
import lsst.afw.detection as afwDetection
import lsst.afw.geom as afwGeom
import lsst.afw.math as afwMath
import lsst.afw.image as afwImage
import lsst.meas.algorithms as measAlg
import lsst.meas.extensions.photometryKron as Kron

try:
    type(verbose)
except NameError:
    verbose = 0
pexLogging.Trace_setVerbosity("meas.photometry.kron", verbose)

try:
    type(display)
except NameError:
    display = False

import lsst.afw.display.ds9 as ds9
import lsst.afw.display.utils as displayUtils

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

class KronPhotometryTestCase(unittest.TestCase):
    """A test case for measuring Kron quantities"""

    def setUp(self):
        pass
        
    def tearDown(self):
        pass

    def makeAndMeasure(self, measureKron, a, b, theta, dx=0.0, dy=0.0, nsigma=6):
        """Make and measure an elliptical Gaussian"""

        width, height = 200, 200
        xcen, ycen = 0.5*width + dx, 0.5*height + dy
        #
        # Make the object
        #
        if a < b:
            a, b = b, a
            theta += 90
        flux = 1e5
        I0 = flux/(2*math.pi*a*b)

        gal = afwImage.ImageF(width, height)

        c, s = math.cos(math.radians(theta)), math.sin(math.radians(theta))
        I, Iuu, Ivv = 0.0, 0.0, 0.0
        for y in range(height):
            for x in range(width):
                dx, dy = x - xcen, y - ycen
                u =  c*dx + s*dy
                v = -s*dx + c*dy
                val = I0*math.exp(-0.5*((u/a)**2 + (v/b)**2))
                if val < 0:
                    val = 0
                gal.set(x, y, val)

                I += val
                Iuu += val*u**2
                Ivv += val*v**2

        Iuu /= I; Ivv /= I

        objImg = afwImage.makeExposure(afwImage.makeMaskedImage(gal))
        objImg.getMaskedImage().getVariance().set(1.0)
        del gal

        kronValues = measureKron(objImg, xcen, ycen, nsigma)

        if display:
            frame = 0
            ds9.mtv(objImg, frame=frame, title="Elliptical")

            ellip = Kron.ellipticalFootprint(afwGeom.makePointI(int(xcen), int(ycen)),
                                             nsigma*a, nsigma*b, math.radians(theta),
                                             afwImage.BBox())

            displayUtils.drawFootprint(ellip, frame=frame)
            ds9.dot("+", xcen, ycen, size=1, ctype=ds9.RED, frame=frame)
            c, s = math.cos(math.radians(theta)), math.sin(math.radians(theta))
            ds9.dot("@:%f,%f,%f" % (nsigma**2*(a**2*c**2 + b**2*s**2),
                                    nsigma**2*(a**2 - b**2)*c*s,
                                    nsigma**2*(a**2*s**2 + b**2*c**2)),
                    xcen, ycen, size=1, ctype=ds9.RED, frame=frame)

        return kronValues

    def measureKron(self, objImg, xcen, ycen, nsigma):
        """Measure Kron quantities using the C++ code"""
        #
        # Now measure some annuli
        #
        mp = measAlg.makeMeasurePhotometry(objImg)
        mp.addAlgorithm("KRON")

        policy = pexPolicy.Policy(pexPolicy.PolicyString(
            """#<?cfg paf policy?>
            KRON: {
                nSigmaForRad: %f
                enabled: true
            }
            """ % (nsigma)))

        mp.configure(policy)
        peak = afwDetection.Peak(xcen, ycen)
        values = mp.measure(peak).find("KRON")
        R_K = values.getParameter()
        flux_K = values.getFlux()
        fluxErr_K = values.getFluxErr()

        return R_K, flux_K, fluxErr_K

    def measureKronInPython(self, objImg, xcen, ycen, nsigma):
        """Measure the Kron quantities of an elliptical Gaussian in python"""
        #
        # Measure moments using SDSS shape algorithm
        #
        ms = measAlg.makeMeasureShape(objImg)
        ms.addAlgorithm("SDSS")

        policy = pexPolicy.Policy(pexPolicy.PolicyString(
            """#<?cfg paf policy?>
            SDSS: {
                enabled: true
            }
            """))

        ms.configure(policy)

        peak = afwDetection.Peak(xcen, ycen)

        values = ms.measure(peak).find("SDSS")
        Mxx = values.getIxx()
        Mxy = values.getIxy()
        Myy = values.getIyy()
        del values
        #
        # Calculate principal axes
        #
        Muu_p_Mvv = Mxx + Myy
        Muu_m_Mvv = math.sqrt((Mxx - Myy)**2 + 4*Mxy**2)
        Muu = 0.5*(Muu_p_Mvv + Muu_m_Mvv)
        Mvv = 0.5*(Muu_p_Mvv - Muu_m_Mvv)
        theta = 0.5*math.atan2(2*Mxy, Mxx - Myy)
        a = math.sqrt(Muu)
        b = math.sqrt(Mvv)
        ab = a/b
        #
        # Get footprint
        #
        ellip = Kron.ellipticalFootprint(afwGeom.makePointI(int(xcen), int(ycen)),
                                         nsigma*a, nsigma*b, theta, afwImage.BBox())
        
        sumI = 0.0
        sumR = 0.38259771140356325/ab*(1 + math.sqrt(2)*math.hypot(math.fmod(xcen, 1), math.fmod(ycen, 1)))*\
               objImg.getMaskedImage().getImage().get(int(xcen), int(ycen))
               
        c, s = math.cos(theta), math.sin(theta)

        for sp in ellip.getSpans():
            y, x0, x1 = sp.getY(), sp.getX0(), sp.getX1()

            for x in range(x0, x1 + 1):
                dx, dy = x - xcen, y - ycen
                u =  c*dx + s*dy
                v = -s*dx + c*dy

                r = math.hypot(u, v*ab)
                val = objImg.getMaskedImage().getImage().get(x, y)

                sumI += val
                sumR += val*r

        R_K = sumR/sumI

        return R_K, 0, 0

    def testEllipticalGaussian(self):
        """Test measuring the Kron quantities of an elliptical Gaussian"""
        #
        # Choose function that does the measuring
        #
        if False:                       # testing only
            measureKron = self.measureKronInPython
        else:
            measureKron = self.measureKron
        #
        # Make and the objects
        #
        ab_vals = (0.5, 1.0, 2.0, 3.0, 4.0, 5.0, )
        #ab_vals = (1.0,)
        for dx in (0.0, 0.5,):
            for dy in (0.0, 0.5,):
                for theta in (20.0, ):
                    for a in ab_vals:
                        for b in ab_vals:
                            if b > a:
                                continue

                            R_K, flux_K, fluxErr_K = self.makeAndMeasure(measureKron, a, b, theta, dx=dx, dy=dy)
                            truth = max(a,b)*math.sqrt(math.pi/2)

                            ID = "a,b %4.1f %4.1f dx,dy = %.1f,%.1f" % (a, b, dx, dy)
                            if not False:
                                print "%s R_K %.3f %.3f %5.1f%%" % (ID, R_K, truth, 100*(R_K/truth - 1))

                            # Set R_K tolerance in pixels
                            if b <= 0.5:
                                if a <= 0.5:
                                    tol = 25
                                elif a <= 1.0:
                                    tol = 15
                                else:
                                    tol = 10
                            elif b <= 1:
                                tol = 2.5
                            else:
                                tol = 0.75

                            if abs(truth - R_K) > 1e-2*tol:
                                self.assertTrue(False,
                                                ("%s  R_Kron: %g v. exact value %g (error %.1f%%)" %
                                                 (ID, R_K, truth, 100*(R_K - truth))))

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def suite():
    """Returns a suite containing all the test cases in this module."""
    tests.init()

    suites = []
    suites += unittest.makeSuite(KronPhotometryTestCase)
    suites += unittest.makeSuite(tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(exit = False):
    """Run the tests"""
    tests.run(suite(), exit)

if __name__ == "__main__":
    run(True)
