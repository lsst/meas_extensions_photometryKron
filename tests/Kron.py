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
import lsst.afw.coord as afwCoord
import lsst.afw.geom as afwGeom
import lsst.afw.geom.ellipses as afwEllipses
import lsst.afw.math as afwMath
import lsst.afw.table as afwTable
import lsst.afw.image as afwImage
import lsst.meas.algorithms as measAlg
import lsst.meas.base as measBase
import lsst.meas.extensions.photometryKron as Kron

try:
    type(verbose)
except NameError:
    verbose = 1
    display = False
    ds9Frame = 0
pexLogging.Trace_setVerbosity("meas.photometry.kron", verbose)

import lsst.afw.display.ds9 as ds9
import lsst.afw.display.utils as displayUtils

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def makeGalaxy(width, height, flux, a, b, theta, dx=0.0, dy=0.0, xy0=None, xcen=None, ycen=None):
    """Make a fake galaxy image"""
    gal = afwImage.ImageF(width, height)
    if xcen is None:
        xcen = 0.5*width + dx
    if ycen is None:
        ycen = 0.5*height + dy
    I0 = flux/(2*math.pi*a*b)

    if xy0 is not None:
        gal.setXY0(xy0)

    c, s = math.cos(math.radians(theta)), math.sin(math.radians(theta))
    I, Iuu, Ivv = 0.0, 0.0, 0.0
    for y in range(height):
        for x in range(width):
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

    exp = afwImage.makeExposure(afwImage.makeMaskedImage(gal))
    exp.getMaskedImage().getVariance().setXY0(exp.getXY0()) # workaround #2577
    exp.getMaskedImage().getVariance().set(1.0)
    exp.setWcs(afwImage.makeWcs(afwCoord.Coord(0.0*afwGeom.degrees, 0.0*afwGeom.degrees),
                                afwGeom.Point2D(0.0, 0.0), 1.0e-4, 0.0, 0.0, 1.0e-4))
    psf = afwDetection.GaussianPsf(15, 15, 3.0)
    exp.setPsf(psf)
    return exp

def makeSourceMeasurementConfig(forced=False, nsigma=6.0, nIterForRadius=1, kfac=2.5):
    """Construct a SourceMeasurementConfig with the requested parameters"""
    if forced:
        msConfig = measBase.ForcedMeasurementConfig()
        msConfig.algorithms.names = ["base_TransformedCentroid", "base_SdssShape", "extensions_photometryKron_KronFlux"]
        msConfig.slots.centroid = "base_TransformedCentroid"
    else:
        msConfig = measBase.SingleFrameMeasurementConfig()
        msConfig.algorithms.names = ["base_SdssCentroid", "base_SdssShape", "extensions_photometryKron_KronFlux"]
        msConfig.slots.centroid = "base_SdssCentroid"
    msConfig.slots.shape = "base_SdssShape"
    msConfig.slots.apFlux = "extensions_photometryKron_KronFlux"
    msConfig.slots.modelFlux = None
    msConfig.slots.psfFlux = None
    msConfig.slots.instFlux = None
    #msConfig.algorithms.names.remove("correctfluxes")
    msConfig.plugins["extensions_photometryKron_KronFlux"].nSigmaForRadius = nsigma
    msConfig.plugins["extensions_photometryKron_KronFlux"].nIterForRadius = nIterForRadius
    msConfig.plugins["extensions_photometryKron_KronFlux"].nRadiusForFlux = kfac
    msConfig.plugins["extensions_photometryKron_KronFlux"].enforceMinimumRadius = False
    return msConfig

def measureFree(exposure, center, msConfig):
    """Unforced measurement"""
    schema = afwTable.SourceTable.makeMinimalSchema()
    task = measBase.SingleFrameMeasurementTask(schema, config=msConfig)
    measCat = afwTable.SourceCatalog(schema)
    source = measCat.addNew()
    ss = afwDetection.FootprintSet(exposure.getMaskedImage(), afwDetection.Threshold(0.1))
    fp = ss.getFootprints()[0]
    source.setFootprint(fp)
    task.run(measCat, exposure)

    return source

# this test should be re-enabled after conversion to meas_base framework on DM-982;
# at present, it uses an old interface that's only on the HSC side
# @unittest.skip("interface unsupported in LSST meas_algorithms")
def measureForced(exposure, source, refWcs, msConfig):
    """Forced measurement"""
    schema = afwTable.SourceTable.makeMinimalSchema()
    task = measBase.ForcedMeasurementTask(schema, config=msConfig)
    forcedCat = task.run(exposure, [source], refWcs).sources
    return forcedCat[0]

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

class KronPhotometryTestCase(tests.TestCase):
    """A test case for measuring Kron quantities"""

    def setUp(self):
        self.flux = 1e5
        self.width, self.height = 200, 200
        self.objImg = None

    def tearDown(self):
        if self.objImg:
            del self.objImg

    def makeAndMeasure(self, measureKron, a, b, theta, dx=0.0, dy=0.0, nsigma=6, kfac=2, nIterForRadius=1,
                       makeImage=True):
        """Make and measure an elliptical Gaussian"""

        xcen, ycen = 0.5*self.width + dx, 0.5*self.height + dy
        #
        # Make the object
        #
        if a < b:
            a, b = b, a
            theta += 90
        I0 = self.flux/(2*math.pi*a*b)

        if makeImage:
            self.objImg = None
        if not self.objImg:
            gal = afwImage.ImageF(self.width, self.height)
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

            self.objImg = afwImage.makeExposure(afwImage.makeMaskedImage(gal))
            self.objImg.getMaskedImage().getVariance().setXY0(self.objImg.getXY0()) # workaround #2577

            self.objImg.getMaskedImage().getVariance().set(1.0)

            if display:
                ds9.mtv(self.objImg, frame=ds9Frame, title="%g %g" % (a, b))

                ds9.dot("+", xcen - self.objImg.getX0(), ycen - self.objImg.getY0(),
                        size=1, ctype=ds9.RED, frame=ds9Frame)
                ds9.pan(xcen - self.objImg.getX0(), ycen - self.objImg.getY0(), frame=ds9Frame)
                c, s = math.cos(math.radians(theta)), math.sin(math.radians(theta))
                # N.b. add 1/12 in quadrature to allow for pixellisation
                ds9.dot("@:%f,%f,%f" % (nsigma**2*((a**2 + 1/12.0)*c**2 + (b**2 + 1/12.0)*s**2),
                                        nsigma**2*(a**2 - b**2)*c*s,
                                        nsigma**2*((a**2 + 1/12.0)*s**2 + (b**2 + 1/12.0)*c**2)),
                        xcen - self.objImg.getX0(), ycen - self.objImg.getY0(),
                        size=1, ctype=ds9.RED, frame=ds9Frame, silent=True)
        #
        # Do the measuring
        #
        FWHM = 5
        ksize = 25                      # size of desired kernel
        self.objImg.setPsf(measAlg.DoubleGaussianPsf(ksize, ksize,
                                                     FWHM/(2*math.sqrt(2*math.log(2))), 1, 0.1))
        return measureKron(self.objImg, xcen, ycen, nsigma, kfac, nIterForRadius)

    def measureKron(self, objImg, xcen, ycen, nsigma, kfac, nIterForRadius):
        """Measure Kron quantities using the C++ code"""
        #
        # Now measure things
        #
        msConfig = measBase.SingleFrameMeasurementConfig()
        msConfig.algorithms.names = ["base_SdssCentroid", "base_SdssShape", "extensions_photometryKron_KronFlux"]
        msConfig.slots.centroid = "base_SdssCentroid"
        msConfig.slots.shape = "base_SdssShape"
        msConfig.slots.apFlux = "extensions_photometryKron_KronFlux"
        msConfig.slots.modelFlux = None
        msConfig.slots.psfFlux = None
        msConfig.slots.instFlux = None
        #msConfig.algorithms.names.remove("correctfluxes")
        msConfig.plugins["extensions_photometryKron_KronFlux"].nSigmaForRadius = nsigma
        msConfig.plugins["extensions_photometryKron_KronFlux"].nIterForRadius = nIterForRadius
        msConfig.plugins["extensions_photometryKron_KronFlux"].nRadiusForFlux = kfac
        msConfig.plugins["extensions_photometryKron_KronFlux"].enforceMinimumRadius = False
        schema = afwTable.SourceTable.makeMinimalSchema()
        task = measBase.SingleFrameMeasurementTask(schema, config=msConfig)
        measCat = afwTable.SourceCatalog(schema)

        source = measCat.addNew()
        ss = afwDetection.FootprintSet(objImg.getMaskedImage(), afwDetection.Threshold(0.1))
        fp = ss.getFootprints()[0]
        source.setFootprint(fp)
        # Then run the default SFM task.  Results not checked
        task.run(measCat, objImg)
        alg = Kron.KronFluxAlgorithm
        R_K = source.get("extensions_photometryKron_KronFlux_radius")
        flux_K = source.get("extensions_photometryKron_KronFlux_flux")
        fluxErr_K = source.get("extensions_photometryKron_KronFlux_fluxSigma")
        flags_K = source.get("extensions_photometryKron_KronFlux_flag_radius")

        if display:
            xc, yc = xcen - objImg.getX0(), ycen - objImg.getY0()
            ds9.dot("x", xc, yc, ctype=ds9.MAGENTA, size=1, frame=ds9Frame)
            displayUtils.drawFootprint(fp, XY0=objImg.getXY0())

            shape = source.getShape()
            if True:                    # nsigma*shape, the radius used to estimate R_K
                shape = shape.clone()
                shape.scale(source.get("extensions_photometryKron_KronFlux_radiusForRadius")/shape.getDeterminantRadius())
                ds9.dot(shape, xc, yc, ctype=ds9.MAGENTA, frame=ds9Frame)
            # Show R_K
            shape = shape.clone()
            for r, ct in [(R_K, ds9.BLUE), (R_K*kfac, ds9.CYAN),]:
                shape.scale(r/shape.getDeterminantRadius())
                ds9.dot(shape, xc, yc, ctype=ct, frame=ds9Frame)

        return R_K, flux_K, fluxErr_K, flags_K

    def measureKronInPython(self, objImg, xcen, ycen, nsigma, kfac, nIterForRadius, makeImage=None):
        """Measure the Kron quantities of an elliptical Gaussian in python

        N.b. only works for XY0 == (0, 0)
        """
        #
        # Measure moments using SDSS shape algorithm
        #
        msConfig = measAlg.SourceMeasurementConfig()
        if False:                       # requires #2546
            msConfig.centroider = None
            msConfig.slots.centroid = None

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

        return R_K, sumI, 0, False

    def xtestEllipticalGaussian(self):
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
        ab_vals = (0.5, 1.0, 2.0, 3.0, 5.0, )
        nIter = 2
        for dx in (0.0, 0.5,):
            for dy in (0.0, 0.5,):
                if measureKron == self.measureKronInPython and dx + dy != 0.0:
                    continue

                for theta in (0.0, 20.0, 45.0, ):
                    for a in ab_vals:
                        for b in ab_vals:
                            if b > a:
                                continue

                            makeImage = True
                            for kfac in (1.5, 2.5,):        # multiple of R_Kron to use for Flux_Kron
                                R_K, flux_K, fluxErr_K, flags_K = \
                                    self.makeAndMeasure(measureKron, a, b, theta, dx=dx, dy=dy, kfac=kfac,
                                                        nIterForRadius=nIter, makeImage=makeImage)
                                makeImage = False
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
                                R_truth = R_truth0*math.sqrt(a*b + 1/12.0*(1 + 0.00286/min(a, b)**3.9))

                                failR = math.isnan(R_K) or flags_K or \
                                    abs(R_truth - R_K) > 1e-2*self.getTolRad(a, b)
                                failFlux =  math.isnan(flux_K) or flags_K or \
                                    abs(flux_K/flux_truth - 1) > 1e-2*self.getTolFlux(a, b, kfac)

                                ID = "a,b,theta %4.1f %4.1f %4.1f  dx,dy = %.1f,%.1f  kfac=%g" % \
                                    (a, b, theta, dx, dy, kfac)
                                if ((failR or failFlux) and verbose) or verbose > 1:
                                    print "%s R_K    %10.3f %10.3f %6.3f pixels (tol %5.3f)%s" % \
                                        (ID, R_K, R_truth, (R_K - R_truth), 1e-2*self.getTolRad(a, b),
                                         " *" if failR else "")
                                    print "%s flux_K %10.3f %10.3f %6.2f%%       (tol %5.3f) %s" % \
                                        (ID, flux_K, flux_truth,
                                         100*(flux_K/flux_truth - 1), self.getTolFlux(a, b, kfac),
                                         " *" if failFlux else "")
                                    if True:
                                        continue # skip tests

                                self.assertFalse(failR, (("%s  R_Kron: %g v. exact value %g " +
                                                          "(error %.3f pixels; limit %.3f)") % \
                                                             (ID, R_K, R_truth, (R_K - R_truth),
                                                              1e-2*self.getTolRad(a, b))))

                                self.assertFalse(failFlux, (("%s  flux_Kron: %g v. exact value %g " +
                                                             "(error %.2f%% limit %.2f%%)") %
                                                            (ID, flux_K, flux_truth, 100*(flux_K/flux_truth-1),
                                                             self.getTolFlux(a, b, kfac))))

    def getTolRad(self, a, b):
        """Return R_K tolerance in hundredths of a pixel"""

        if b <= 0.5:
            if a <= 0.5:
                tol = 35
            elif a <= 2:
                tol = 350
            else:
                tol = 25*a              # i.e. 0.25*a
        elif b <= 1:
            tol = 7.0
        else:
            tol = 1.0

        return tol

    def getTolFlux(self, a, b, kfac):
        """Return Flux_K tolerance in percent"""
        if b <= 0.5:
            if a <= 0.5:
                if kfac > 2:
                    tol = 5.0
                else:
                    tol = 10.0
            elif a <= 1.0:
                if kfac <= 1.5:
                    tol = 10.0
                else:
                    tol = 4.0
            else:
                if kfac > 2:
                    tol = 3.0
                elif kfac > 1.5:
                    tol = 5.0
                else:
                    tol = 10.0
        elif b <= 1:
            if a <= 1:
                tol = 2.0
            else:
                if kfac > 2:
                    tol = 0.25
                elif kfac > 1.5:
                    tol = 0.5
                else:
                    tol = 1.27
        elif b <= 2:
            if kfac > 1.5:
                tol = 0.1
            else:
                tol = 0.5
        else:
            tol = 0.30

        return tol

    def testForced(self):
        """Check that forced photometry works in the presence of rotations and translations"""
        kfac = 2.5
        warper = afwMath.Warper("lanczos4")
        a = 13
        for axisRatio in (0.25, 1.0):
            b = a*axisRatio
            for theta in (0, 30, 45):
                width, height = 256, 256
                center = afwGeom.Point2D(0.5*width, 0.5*height)
                original = makeGalaxy(width, height, 1000.0, a, b, theta)
                msConfig = makeSourceMeasurementConfig(kfac=kfac)
                source = measureFree(original, center, msConfig)
                if source.get("extensions_photometryKron_KronFlux_flag"):
                    continue
                angleList = [45, 90,]
                scaleList = [1.0, 0.5]
                offsetList = [(1.23, 4.56), (12.3, 45.6)]

                for angle in angleList:
                 for scale in scaleList:
                  for offset in offsetList:
#angle, scale, offset in itertools.product(angleList, scaleList, offsetList):
                    cosAngle = math.cos(math.radians(angle))
                    sinAngle = math.sin(math.radians(angle))
                    dx, dy = offset
                    pixelScale = original.getWcs().pixelScale().asDegrees()*scale
                    wcs = afwImage.makeWcs(afwCoord.Coord(0.0*afwGeom.degrees, 0.0*afwGeom.degrees),
                                           afwGeom.Point2D(dx, dy), pixelScale*cosAngle,
                                           pixelScale*sinAngle, -pixelScale*sinAngle, pixelScale*cosAngle)

                    warped = warper.warpExposure(wcs, original)
                    msConfig = makeSourceMeasurementConfig(kfac=kfac, forced=True)
                    forced = measureForced(warped, source, original.getWcs(), msConfig)

                    if display:
                        ds9.mtv(original, frame=1)
                        shape = source.getShape().clone()
                        xc, yc = source.getCentroid()
                        radius = source.get("extensions_photometryKron_KronFlux_radius")
                        for r, ct in [(radius, ds9.BLUE), (radius*kfac, ds9.CYAN),]:
                            shape.scale(r/shape.getDeterminantRadius())
                            ds9.dot(shape, xc, yc, ctype=ct, frame=1)
                        ds9.mtv(warped, frame=2)
                        transform = (wcs.linearizeSkyToPixel(forced.getCoord())*
                                     original.getWcs().linearizePixelToSky(source.getCoord()))
                        shape = shape.transform(transform.getLinear())
                        radius = forced.get("extensions_photometryKron_KronFlux_radius")
                        xc, yc = wcs.skyToPixel(source.getCoord()) - afwGeom.Extent2D(warped.getXY0())
                        for r, ct in [(radius, ds9.BLUE), (radius*kfac, ds9.CYAN),]:
                            shape.scale(r/shape.getDeterminantRadius())
                            ds9.dot(shape, xc, yc, ctype=ct, frame=2)

                    try:

                        #self.assertClose(source.get("extensions_photometryKron_KronFlux_flux"), forced.get("extensions_photometryKron_KronFlux_flux"), rtol=2.0e-4, atol=None)
                        #self.assertClose(source.get("extensions_photometryKron_KronFlux_radius"), scale*forced.get("extensions_photometryKron_KronFlux_radius"), rtol=None, atol=1.0e-12)
                        self.assertEqual(source.get("extensions_photometryKron_KronFlux_flag"), forced.get("extensions_photometryKron_KronFlux_flag"))
                    except:
                        print ("Failed:", angle, scale, offset,)
                             #  [(source.get(f), forced.get(f)) for f in
                             #   ("flux.kron", "flux.kron.radius", "flux.kron.flags")])
                        raise


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
