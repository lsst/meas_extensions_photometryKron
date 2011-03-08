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
import lsst.afw.math as afwMath
import lsst.afw.image as afwImage
import lsst.meas.algorithms as measAlg
import lsst.meas.extensions.photometryKron

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

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

class KronPhotometryTestCase(unittest.TestCase):
    """A test case for measuring Kron quantities"""

    def setUp(self):
        pass
        
    def tearDown(self):
        pass

    def testEllipticalGaussian(self):
        """Test measuring the properties of an elliptical Gaussian"""

        width, height = 200, 200
        xcen, ycen = 0.5*width, 0.5*height
        #
        # Make the object
        #
        gal = afwImage.ImageF(width, height)
        a, b, theta = float(10), float(5), 20
        flux = 1e4
        I0 = flux/(2*math.pi*a*b)

        c, s = math.cos(math.radians(theta)), math.sin(math.radians(theta))
        for y in range(height):
            for x in range(width):
                dx, dy = x - xcen, y - ycen
                u =  c*dx + s*dy
                v = -s*dx + c*dy
                val = I0*math.exp(-0.5*((u/a)**2 + (v/b)**2))
                if val < 0:
                    val = 0
                gal.set(x, y, val)

        objImg = afwImage.makeExposure(afwImage.makeMaskedImage(gal))
        objImg.getMaskedImage().getVariance().set(1.0)
        del gal

        if display:
            frame = 0
            ds9.mtv(objImg, frame=frame, title="Elliptical")

        self.assertAlmostEqual(1.0, afwMath.makeStatistics(objImg.getMaskedImage().getImage(),
                                                           afwMath.SUM).getValue()/flux)
        #
        # Now measure some annuli
        #
        mp = measAlg.makeMeasurePhotometry(objImg)
        mp.addAlgorithm("KRON")
    
        policy = pexPolicy.Policy(pexPolicy.PolicyString(
            """#<?cfg paf policy?>
            KRON: {
                enabled: true
            }
            """))
        
        mp.configure(policy)

        peak = afwDetection.Peak(xcen, ycen)
        photom = mp.measure(peak)

        values = photom.find("KRON")
        print photom.get("KRON", "RADIUS")
        print values.getRadius(0), values.getFlux(), values.getFluxErr()

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
