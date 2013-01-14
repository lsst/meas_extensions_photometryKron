#!/usr/bin/env python
"""
Tests for measuring things

Run with:
   python Kron.py
or
   python
   >>> import Kron; Kron.run()
"""

import math, os, sys, unittest
import numpy as np
import lsst.utils.tests as tests
import lsst.pex.exceptions as pexExceptions
import lsst.pex.logging as pexLogging
import lsst.pex.policy as pexPolicy
import lsst.afw.detection as afwDetection
import lsst.afw.geom as afwGeom
import lsst.afw.geom.ellipses as afwEllipses
import lsst.afw.math as afwMath
import lsst.afw.table as afwTable
import lsst.afw.image as afwImage
import lsst.meas.algorithms as measAlg
import lsst.meas.extensions.photometryKron as Kron

try:
    type(verbose)
except NameError:
    verbose = 0
pexLogging.Trace_setVerbosity("meas.photometry.kron", verbose)

try:
    display
except NameError:
    ds9Frame = 0
    display = False

import lsst.afw.display.ds9 as ds9
import lsst.afw.display.utils as displayUtils

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

class KronPhotometryTestCase(unittest.TestCase):
    """A test case for measuring Kron quantities"""

    def setUp(self):
        self.flux = 1e5
        self.width, self.height = 200, 200
        
    def tearDown(self):
        pass

    def makeAndMeasure(self, measureKron, a, b, theta, dx=0.0, dy=0.0, nsigma=6, kfac=2):
        """Make and measure an elliptical Gaussian"""

        xcen, ycen = 0.5*self.width + dx, 0.5*self.height + dy
        #
        # Make the object
        #
        if a < b:
            a, b = b, a
            theta += 90
        I0 = self.flux/(2*math.pi*a*b)

        gal = afwImage.ImageF(self.width, self.height)
        if True:
            global warnedXY0
            try: warnedXY0
            except:            
                print "Not setting XY0"
                warnedXY0 = True
        else:
            gal.setXY0(10, 10)

        c, s = math.cos(math.radians(theta)), math.sin(math.radians(theta))
        I, Iuu, Ivv = 0.0, 0.0, 0.0
        for y in range(self.height):
            for x in range(self.width):
                dx, dy = x + gal.getX0() - xcen, y + gal.getY0() - ycen
                if math.hypot(dx, dy) < 10.5:
                    nsample = float(5)
                    subZ = np.linspace(-0.5*(1 - 1/nsample), 0.5*(1 - 1/nsample), nsample)
                else:
                    nsample = 1
                    subZ = [0.0]

                val = 0
                for sx in subZ:
                    for sy in subZ:
                        u =  c*(dx + sx) + s*(dy + sy)
                        v = -s*(dx + sx) + c*(dy + sy)
                        val += I0*math.exp(-0.5*((u/a)**2 + (v/b)**2))

                if val < 0:
                    val = 0
                gal.set(x, y, val/nsample**2)

                I += val
                Iuu += val*u**2
                Ivv += val*v**2

        Iuu /= I; Ivv /= I

        objImg = afwImage.makeExposure(afwImage.makeMaskedImage(gal))
        objImg.getMaskedImage().getVariance().set(1.0)

        FWHM = 5
        ksize = 25                      # size of desired kernel
        objImg.setPsf(afwDetection.createPsf("DoubleGaussian", ksize, ksize,
                                             FWHM/(2*math.sqrt(2*math.log(2))), 1, 0.1))

        if display:
            ds9.mtv(objImg, frame=ds9Frame, title="Elliptical")

            ellipse = afwEllipses.Ellipse(afwEllipses.Axes(nsigma*a, nsigma*b, theta),
                                          afwGeom.Point2D(xcen - objImg.getX0(), ycen - objImg.getY0()))
            fpEllipse = afwDetection.Footprint(ellipse)

            displayUtils.drawFootprint(fpEllipse, frame=ds9Frame)
            ds9.dot("+", xcen - gal.getX0(), ycen - gal.getY0(), size=1, ctype=ds9.RED, frame=ds9Frame)
            ds9.pan(xcen - gal.getX0(), ycen - gal.getY0(), frame=ds9Frame)
            c, s = math.cos(math.radians(theta)), math.sin(math.radians(theta))
            ds9.dot("@:%f,%f,%f" % (nsigma**2*(a**2*c**2 + b**2*s**2),
                                    nsigma**2*(a**2 - b**2)*c*s,
                                    nsigma**2*(a**2*s**2 + b**2*c**2)),
                    xcen - gal.getX0(), ycen - gal.getY0(), size=1, ctype=ds9.RED, frame=ds9Frame)
        #
        # Do the measuring
        #
        return measureKron(objImg, xcen, ycen, nsigma, kfac)

    def measureKron(self, objImg, xcen, ycen, nsigma, kfac):
        """Measure Kron quantities using the C++ code"""
        #
        # Now measure things
        #
        msConfig = measAlg.SourceMeasurementConfig()
        msConfig.algorithms.names.add("flux.kron")
        msConfig.algorithms["flux.kron"].nSigmaForRadius = nsigma
        msConfig.algorithms["flux.kron"].nRadiusForFlux = kfac
        schema = afwTable.SourceTable.makeMinimalSchema()
        ms = msConfig.makeMeasureSources(schema)
        table = afwTable.SourceTable.make(schema)
        msConfig.slots.setupTable(table)
        source = table.makeRecord()
        fp = afwDetection.Footprint(objImg.getBBox())
        source.setFootprint(fp)
        center = afwGeom.Point2D(xcen, ycen)
        ms.apply(source, objImg, center)

        R_K = source.get("flux.kron.radius")
        flux_K = source.get("flux.kron")
        fluxErr_K = source.get("flux.kron.err")

        if display:
            shape = source.getShape()
            if True:                    # nsigma*shape, the radius used to estimate R_K
                shape = shape.clone()
                shape.scale(nsigma)
                ds9.dot("@:%f,%f,%f" % (shape.getIxx(), shape.getIxy(), shape.getIyy()), 
                        *(source.getCentroid() - afwGeom.Extent2D(objImg.getXY0())),
                        ctype=ds9.MAGENTA, frame=ds9Frame)
            # Show R_K
            shape = shape.clone(); shape.scale(R_K/shape.getTraceRadius())
            ds9.dot("@:%f,%f,%f" % (shape.getIxx(), shape.getIxy(), shape.getIyy()), 
                    *(source.getCentroid() - afwGeom.Extent2D(objImg.getXY0())),
                    ctype=ds9.BLUE, frame=ds9Frame)

        return R_K, flux_K, fluxErr_K

    def measureKronInPython(self, objImg, xcen, ycen, nsigma, kfac):
        """Measure the Kron quantities of an elliptical Gaussian in python

        N.b. only works for XY0 == (0, 0)
        """
        #
        # Measure moments using SDSS shape algorithm
        #
        #
        # Now measure things
        #
        msConfig = measAlg.SourceMeasurementConfig()
        schema = afwTable.SourceTable.makeMinimalSchema()
        ms = msConfig.makeMeasureSources(schema)
        table = afwTable.SourceTable.make(schema)
        msConfig.slots.setupTable(table)
        source = table.makeRecord()
        fp = afwDetection.Footprint(objImg.getBBox())
        source.setFootprint(fp)
        center = afwGeom.Point2D(xcen, ycen)
        ms.apply(source, objImg, center)

        Mxx = source.getIxx()
        Mxy = source.getIxy()
        Myy = source.getIyy()
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
        ellipse = afwEllipses.Ellipse(afwEllipses.Axes(nsigma*a, nsigma*b, theta),
                                      afwGeom.Point2D(xcen - objImg.getX0(), ycen - objImg.getY0()))
        fpEllipse = afwDetection.Footprint(ellipse)
        
        sumI = 0.0
        sumR = 0.38259771140356325/ab*(1 + math.sqrt(2)*math.hypot(math.fmod(xcen, 1), math.fmod(ycen, 1)))*\
               objImg.getMaskedImage().getImage().get(int(xcen), int(ycen))
               
        gal = objImg.getMaskedImage().getImage()

        c, s = math.cos(theta), math.sin(theta)
        for sp in fpEllipse.getSpans():
            y, x0, x1 = sp.getY(), sp.getX0(), sp.getX1()

            for x in range(x0, x1 + 1):
                dx, dy = x - xcen, y - ycen
                u =  c*dx + s*dy
                v = -s*dx + c*dy

                r = math.hypot(u, v*ab)
                try:
                    val = gal.get(x, y)
                except:
                    continue

                sumI += val
                sumR += val*r

        R_K = sumR/sumI

        sumI = 0.0
        for y in range(self.height):
            for x in range(self.width):
                dx, dy = x - xcen, y - ycen
                u =  c*dx + s*dy
                v = -s*dx + c*dy
                if math.hypot(u/a, v/b) < kfac:
                    sumI += gal.get(x, y)

        return R_K, sumI, 0

    def testEllipticalGaussian(self):
        """Test measuring the Kron quantities of an elliptical Gaussian"""
        #
        # Choose function that does the measuring
        #
        if False:                       # testing only; requires XY0 == (0, 0)
            measureKron = self.measureKronInPython
        else:
            measureKron = self.measureKron
        #
        # Make and the objects
        #
        kfac = 2.5                      # multiple of R_Kron to use for Flux_Kron
        kfac = 1.5

        ab_vals = (0.5, 1.0, 2.0, 3.0, 4.0, 5.0, )
        for dx in (0.0, 0.5,):
            for dy in (0.0, 0.5,):
                if measureKron == self.measureKronInPython and dx + dy != 0.0:
                    continue
                
                for theta in (20.0, ):
                    for a in ab_vals:
                        for b in ab_vals:
                            if b > a:
                                continue

                            R_K, flux_K, fluxErr_K = self.makeAndMeasure(measureKron, a, b, theta,
                                                                         dx=dx, dy=dy, kfac=kfac)
                            #
                            # We'll have to correct for the pixelisation as we sum over the central
                            # few pixels when making models, mostly do deal with b ~ 0.5 models.
                            #
                            # See Section 5 of
                            #   http://www.astro.princeton.edu/~rhl/photomisc/aperture.pdf
                            # for the source of 0.00286 etc.
                            #
                            R_truth0 = math.sqrt(math.pi/2)
                            R_truth = R_truth0*math.sqrt(1 + 0.8*1/(12.0*a*b))

                            flux_truth = self.flux*(1 - math.exp(-0.5*(kfac*R_truth)**2))
                            R_truth = R_truth0*math.sqrt(max(a,b)**2 + 1/12.0*(1 + 0.00286/min(a, b)**3.9))
                            ID = "a,b %4.1f %4.1f dx,dy = %.1f,%.1f" % (a, b, dx, dy)
                            if False:
                                print "%s R_K %.3f %.3f %5.2f pixels" % (ID, R_K, R_truth, (R_K - R_truth))
                            if False:
                                print "%s flux_K %.3f %.3f %5.1f%%" % (ID, flux_K, flux_truth,
                                                                       100*(flux_K/flux_truth - 1))

                            if math.isnan(R_K) or \
                                    abs(R_truth - R_K) > 1e-2*self.getTolRad(a, b):
                                self.assertTrue(False,
                                                ("%s  R_Kron: %g v. exact value %g (error %.2f pixels)" %
                                                 (ID, R_K, R_truth, (R_K - R_truth))))
                            if math.isnan(flux_K) or \
                                    abs(flux_truth/flux_K - 1) > 1e-2*self.getTolFlux(a, b, kfac):
                                self.assertTrue(False,
                                                ("%s  flux_Kron: %g v. exact value %g (error %.1f%%)" %
                                                 (ID, flux_K, flux_truth, 100*(flux_K/flux_truth - 1))))

    def getTolRad(self, a, b):
        """Return R_K tolerance in hundredths of a pixel"""

        if b <= 0.5:
            if a <= 0.5:
                tol = 25
            elif a <= 1:
                tol = 14
            elif a <= 2:
                tol = 10
            else:
                tol = 6
        elif b <= 1:
            tol = 2.5
        else:
            tol = 0.75

        return tol

    def getTolFlux(self, a, b, kfac):
        """Return Flux_K tolerance in percent"""
        if b <= 0.5:
            if a <= 0.5:
                if kfac > 2:
                    tol = 5.0
                else:
                    tol = 16.0
            else:
                if kfac > 2:
                    tol = 3.0
                elif kfac > 1.5:
                    tol = 5.0
                else:
                    tol = 7.0
        elif b <= 1:
            if kfac > 2:
                tol = 0.25
            elif kfac > 1.5:
                tol = 0.5
            else:
                tol = 1.2
        elif b <= 2:
            if kfac > 1.5:
                tol = 0.1
            else:
                tol = 0.6
        else:
            tol = 0.1

        return tol

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
