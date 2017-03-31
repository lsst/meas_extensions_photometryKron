#!/usr/bin/env python
#
# LSST Data Management System
#
# Copyright 2008-2016  AURA/LSST.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <https://www.lsstcorp.org/LegalNotices/>.
#
from __future__ import print_function
from builtins import range
import math
import unittest
import sys

import numpy as np
import itertools
import lsst.utils.tests
import lsst.afw.detection as afwDetection
import lsst.afw.geom as afwGeom
import lsst.afw.geom.ellipses as afwEllipses
import lsst.afw.coord as afwCoord
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.afw.table as afwTable
import lsst.meas.algorithms as measAlg
import lsst.meas.base as measBase
# importing this package registers essential code
import lsst.meas.extensions.photometryKron
from lsst.daf.base import PropertyList

try:
    type(verbose)
except NameError:
    display = False
    ds9Frame = 0

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
                    u = c*(dx + sx) + s*(dy + sy)
                    v = -s*(dx + sx) + c*(dy + sy)
                    val += I0*math.exp(-0.5*((u/a)**2 + (v/b)**2))

            if val < 0:
                val = 0
            gal.set(x, y, val/nsample**2)

            I += val
            Iuu += val*u**2
            Ivv += val*v**2

    Iuu /= I
    Ivv /= I

    exp = afwImage.makeExposure(afwImage.makeMaskedImage(gal))
    exp.getMaskedImage().getVariance().set(1.0)
    exp.setWcs(afwImage.makeWcs(afwCoord.Coord(0.0*afwGeom.degrees, 0.0*afwGeom.degrees),
                                afwGeom.Point2D(0.0, 0.0), 1.0e-4, 0.0, 0.0, 1.0e-4))
    # add a dummy Psf.  The new SdssCentroid needs one
    exp.setPsf(afwDetection.GaussianPsf(11, 11, 0.01))
    return exp


def makeMeasurementConfig(forced=False, nsigma=6.0, nIterForRadius=1, kfac=2.5):
    """Construct a (SingleFrame|Forced)MeasurementConfig with the requested parameters"""
    if forced:
        msConfig = measBase.ForcedMeasurementConfig()
        msConfig.algorithms.names = ["base_TransformedCentroid", "base_TransformedShape",
                                     "ext_photometryKron_KronFlux"]
        msConfig.slots.centroid = "base_TransformedCentroid"
        msConfig.slots.shape = "base_TransformedShape"
        msConfig.copyColumns = {"id": "objectId", "parent": "parentObjectId"}
    else:
        msConfig = measBase.SingleFrameMeasurementConfig()
        msConfig.algorithms.names = ["base_SdssCentroid", "base_SdssShape",
                                     "ext_photometryKron_KronFlux", "base_SkyCoord"]
        msConfig.slots.centroid = "base_SdssCentroid"
        msConfig.slots.shape = "base_SdssShape"
    msConfig.slots.apFlux = "ext_photometryKron_KronFlux"
    msConfig.slots.modelFlux = None
    msConfig.slots.psfFlux = None
    msConfig.slots.instFlux = None
    msConfig.slots.calibFlux = None
    # msConfig.algorithms.names.remove("correctfluxes")
    msConfig.plugins["ext_photometryKron_KronFlux"].nSigmaForRadius = nsigma
    msConfig.plugins["ext_photometryKron_KronFlux"].nIterForRadius = nIterForRadius
    msConfig.plugins["ext_photometryKron_KronFlux"].nRadiusForFlux = kfac
    msConfig.plugins["ext_photometryKron_KronFlux"].enforceMinimumRadius = False
    return msConfig


def measureFree(exposure, center, msConfig):
    """Unforced measurement"""
    schema = afwTable.SourceTable.makeMinimalSchema()
    algMeta = PropertyList()
    task = measBase.SingleFrameMeasurementTask(schema, config=msConfig, algMetadata=algMeta)
    measCat = afwTable.SourceCatalog(schema)
    source = measCat.addNew()
    source.getTable().setMetadata(algMeta)
    ss = afwDetection.FootprintSet(exposure.getMaskedImage(), afwDetection.Threshold(0.1))
    fp = ss.getFootprints()[0]
    source.setFootprint(fp)
    task.run(measCat, exposure)
    return source


def measureForced(exposure, source, refWcs, msConfig):
    """Forced measurement"""
    refCat = afwTable.SourceCatalog(source.table)
    refCat.append(source)
    schema = afwTable.SourceTable.makeMinimalSchema()
    algMeta = PropertyList()
    source.getTable().setMetadata(algMeta)
    task = measBase.ForcedMeasurementTask(schema, config=msConfig, algMetadata=algMeta)
    measCat = task.generateMeasCat(exposure, refCat, refWcs)
    task.attachTransformedFootprints(measCat, refCat, exposure, refWcs)
    task.run(measCat, exposure, refCat, refWcs)
    return measCat[0]


class KronPhotometryTestCase(lsst.utils.tests.TestCase):
    """A test case for measuring Kron quantities"""

    def setUp(self):
        self.flux = 1e5
        self.width, self.height = 200, 200
        self.objImg = None

    def tearDown(self):
        if self.objImg:
            del self.objImg

    def makeAndMeasure(self, measureKron, a, b, theta, dx=0.0, dy=0.0, nsigma=6, kfac=2, nIterForRadius=1,
                       xcen=None, ycen=None,
                       makeImage=True):
        """Make and measure an elliptical Gaussian"""

        if xcen is None:
            xcen = 0.5*self.width + dx
        if ycen is None:
            ycen = 0.5*self.height + dy
        #
        # Make the object
        #
        if a < b:
            a, b = b, a
            theta += 90

        if self.objImg is None:
            makeImage = True
        if makeImage:
            self.objImg = makeGalaxy(self.width, self.height, self.flux, a, b, theta, dx, dy,
                                     afwGeom.Point2I(10, 10), xcen=xcen, ycen=ycen)

            if display:
                ds9.mtv(self.objImg, frame=ds9Frame, title="%g %g" % (a, b))
        if display:
            if not makeImage:
                ds9.erase(frame=ds9Frame)

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
        center = afwGeom.Point2D(xcen, ycen)
        msConfig = makeMeasurementConfig(False, nsigma, nIterForRadius, kfac)
        source = measureFree(objImg, center, msConfig)
        algMeta = source.getTable().getMetadata()
        self.assertTrue(algMeta.exists('ext_photometryKron_KronFlux_nRadiusForFlux'))

        R_K = source.get("ext_photometryKron_KronFlux_radius")
        flux_K = source.get("ext_photometryKron_KronFlux_flux")
        fluxErr_K = source.get("ext_photometryKron_KronFlux_fluxSigma")
        flags_K = source.get("ext_photometryKron_KronFlux_flag")
        if not flags_K:
            # Forced measurement on the same image should produce exactly the same result
            msConfig = makeMeasurementConfig(True, nsigma, nIterForRadius, kfac)
            forced = measureForced(objImg, source, objImg.getWcs(), msConfig)
            algMeta = source.getTable().getMetadata()
            self.assertTrue(algMeta.exists('ext_photometryKron_KronFlux_nRadiusForFlux'))
            for field in (
                "ext_photometryKron_KronFlux_flux",
                "ext_photometryKron_KronFlux_fluxSigma",
                "ext_photometryKron_KronFlux_radius",
                "ext_photometryKron_KronFlux_flag"
            ):
                try:
                    if np.isnan(source.get(field)):
                        self.assertTrue(np.isnan(forced.get(field)))
                    else:
                        self.assertFloatsAlmostEqual(source.get(
                            field), forced.get(field), rtol=1.0e-6, atol=None)
                except AssertionError:
                    print("Failed:", field, source.get(field), forced.get(field))
                    raise

        if display:
            xc, yc = xcen - objImg.getX0(), ycen - objImg.getY0()
            ds9.dot("x", xc, yc, ctype=ds9.MAGENTA, size=1, frame=ds9Frame)
            displayUtils.drawFootprint(source.getFootprint(), XY0=objImg.getXY0())

            shape = source.getShape()
            if True:                    # nsigma*shape, the radius used to estimate R_K
                shape = shape.clone()
                shape.scale(source.get("ext_photometryKron_KronFlux_radiusForRadius") /
                            shape.getDeterminantRadius())
                ds9.dot(shape, xc, yc, ctype=ds9.MAGENTA, frame=ds9Frame)
            # Show R_K
            shape = shape.clone()
            for r, ct in [(R_K, ds9.BLUE), (R_K*kfac, ds9.CYAN), ]:
                shape.scale(r/shape.getDeterminantRadius())
                ds9.dot(shape, xc, yc, ctype=ct, frame=ds9Frame)

        return R_K, flux_K, fluxErr_K, flags_K, \
            source.get("ext_photometryKron_KronFlux_flag_bad_radius"), \
            source.get("ext_photometryKron_KronFlux_flag_small_radius")

    def measureKronInPython(self, objImg, xcen, ycen, nsigma, kfac, nIterForRadius, makeImage=None):
        """Measure the Kron quantities of an elliptical Gaussian in python

        N.b. only works for XY0 == (0, 0)
        """
        #
        # Measure moments using SDSS shape algorithm
        #
        # Note: this code was converted to the new meas_framework, but is not exercised.
        msConfig = makeMeasurementConfig(False, nsigma, nIterForRadius, kfac)
        center = afwGeom.Point2D(xcen, ycen)
        source = self.measureFree(objImg, center, msConfig)
        algMeta = source.getTable().getMetadata()
        self.assertTrue(algMeta.exists('ext_photometryKron_KronFlux_nRadiusForFlux'))

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
        sumR = 0.38259771140356325/ab*(1 + math.sqrt(2)*math.hypot(math.fmod(xcen, 1), math.fmod(ycen, 1))) *\
            objImg.getMaskedImage().getImage().get(int(xcen), int(ycen))

        gal = objImg.getMaskedImage().getImage()

        c, s = math.cos(theta), math.sin(theta)
        for sp in fpEllipse.getSpans():
            y, x0, x1 = sp.getY(), sp.getX0(), sp.getX1()

            for x in range(x0, x1 + 1):
                dx, dy = x - xcen, y - ycen
                u = c*dx + s*dy
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
                u = c*dx + s*dy
                v = -s*dx + c*dy
                if math.hypot(u/a, v/b) < kfac:
                    sumI += gal.get(x, y)

        return R_K, sumI, 0, False, False, False

    def testEllipticalGaussian(self):
        """Test measuring the Kron quantities of an elliptical Gaussian"""

        ignoreTestFailures = False  # if True, keep going after test failures but always generate a failure
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
                                R_K, flux_K, fluxErr_K, flags_K, flags_radius, flags_small_radius = \
                                    self.makeAndMeasure(measureKron, a, b, theta, dx=dx, dy=dy, kfac=kfac,
                                                        nIterForRadius=nIter, makeImage=makeImage)
                                makeImage = False

                                self.assertFalse(flags_radius)
                                self.assertFalse(flags_small_radius)
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
                                    print("%s R_K    %10.3f %10.3f %6.3f pixels (tol %5.3f)%s" % \
                                        (ID, R_K, R_truth, (R_K - R_truth), 1e-2*self.getTolRad(a, b),
                                         " *" if failR else ""))
                                    print("%s flux_K %10.3f %10.3f %6.2f%%       (tol %5.3f) %s" % \
                                        (ID, flux_K, flux_truth,
                                         100*(flux_K/flux_truth - 1), self.getTolFlux(a, b, kfac),
                                         " *" if failFlux else ""))

                                if ignoreTestFailures:
                                    continue

                                self.assertFalse(failR, (("%s  R_Kron: %g v. exact value %g " +
                                                          "(error %.3f pixels; limit %.3f)") %
                                                         (ID, R_K, R_truth, (R_K - R_truth),
                                                          1e-2*self.getTolRad(a, b))))

                                self.assertFalse(failFlux,
                                                 (("%s  flux_Kron: %g v. exact value %g " +
                                                   "(error %.2f%% limit %.2f%%)") %
                                                  (ID, flux_K, flux_truth, 100*(flux_K/flux_truth-1),
                                                   self.getTolFlux(a, b, kfac))))

        self.assertFalse(ignoreTestFailures, "You are ignoring possible test failures")

    def testBadFlags(self):
        a, b, theta = 5, 1, 20.0
        #
        # Check that we fail if too close to the edge of the image
        #
        for cen in (10, 20, 38, 40, 50, self.width - 10):
            makeImage = True
            for kfac in (2.5, 10,):
                R_K, flux_K, fluxErr_K, flags_K, flags_radius, flags_small_radius = \
                    self.makeAndMeasure(self.measureKron, a, b, theta, makeImage=makeImage,
                                        xcen=cen, ycen=cen, kfac=kfac)
                makeImage = False

                msg = "KronFlux_flag_bad_radius: cen = (%g, %g), kfac = %g" % (cen, cen, kfac)

                if kfac == 2.5 and (cen <= 20 or cen > self.width - 20):
                    self.assertTrue(flags_K, msg)
                elif kfac == 10:
                    self.assertTrue(flags_K, msg)
                else:
                    self.assertFalse(flags_K, msg)

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
                msConfig = makeMeasurementConfig(forced=False, kfac=kfac)
                source = measureFree(original, center, msConfig)
                algMeta = source.getTable().getMetadata()
                self.assertTrue(algMeta.exists('ext_photometryKron_KronFlux_nRadiusForFlux'))
                if source.get("ext_photometryKron_KronFlux_flag"):
                    continue

                angleList = [45, 90, ]
                scaleList = [1.0, 0.5]
                offsetList = [(1.23, 4.56), (12.3, 45.6)]

                for angle, scale, offset in itertools.product(angleList, scaleList, offsetList):
                    cosAngle = math.cos(math.radians(angle))
                    sinAngle = math.sin(math.radians(angle))
                    dx, dy = offset
                    pixelScale = original.getWcs().pixelScale().asDegrees()*scale
                    wcs = afwImage.makeWcs(afwCoord.Coord(0.0*afwGeom.degrees, 0.0*afwGeom.degrees),
                                           afwGeom.Point2D(dx, dy), pixelScale*cosAngle,
                                           pixelScale*sinAngle, -pixelScale*sinAngle, pixelScale*cosAngle)

                    warped = warper.warpExposure(wcs, original)
                    # add a Psf if there is none.  The new SdssCentroid needs a Psf.
                    if warped.getPsf() == None:
                        warped.setPsf(afwDetection.GaussianPsf(11, 11, 0.01))
                    msConfig = makeMeasurementConfig(kfac=kfac, forced=True)
                    forced = measureForced(warped, source, original.getWcs(), msConfig)
                    algMeta = source.getTable().getMetadata()
                    self.assertTrue(algMeta.exists('ext_photometryKron_KronFlux_nRadiusForFlux'))

                    if display:
                        ds9.mtv(original, frame=1)
                        shape = source.getShape().clone()
                        xc, yc = source.getCentroid()
                        radius = source.get("ext_photometryKron_KronFlux_radius")
                        for r, ct in [(radius, ds9.BLUE), (radius*kfac, ds9.CYAN), ]:
                            shape.scale(r/shape.getDeterminantRadius())
                            ds9.dot(shape, xc, yc, ctype=ct, frame=1)
                        ds9.mtv(warped, frame=2)
                        transform = (wcs.linearizeSkyToPixel(source.getCoord()) *
                                     original.getWcs().linearizePixelToSky(source.getCoord()))
                        shape = shape.transform(transform.getLinear())
                        radius = source.get("ext_photometryKron_KronFlux_radius")
                        xc, yc = wcs.skyToPixel(source.getCoord()) - afwGeom.Extent2D(warped.getXY0())
                        for r, ct in [(radius, ds9.BLUE), (radius*kfac, ds9.CYAN), ]:
                            shape.scale(r/shape.getDeterminantRadius())
                            ds9.dot(shape, xc, yc, ctype=ct, frame=2)
                    try:
                        self.assertFloatsAlmostEqual(source.get("ext_photometryKron_KronFlux_flux"),
                                                     forced.get("ext_photometryKron_KronFlux_flux"),
                                                     rtol=1.0e-3
                                                     )
                        self.assertFloatsAlmostEqual(source.get("ext_photometryKron_KronFlux_radius"),
                                                     scale*forced.get("ext_photometryKron_KronFlux_radius"),
                                                     rtol=1.0e-3
                                                     )
                        self.assertEqual(source.get("ext_photometryKron_KronFlux_flag"),
                                         forced.get("ext_photometryKron_KronFlux_flag")
                                         )
                    except:
                        print(("Failed:", angle, scale, offset,
                              [(source.get(f), forced.get(f)) for f in
                               ("ext_photometryKron_KronFlux_flux",
                                "ext_photometryKron_KronFlux_radius",
                                "ext_photometryKron_KronFlux_flag"
                                )
                               ]))
                        raise


#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("--display", default=False, action="store_true", help="Activate display?")
    parser.add_argument("--verbose", type=int, default=0, help="Verbosity level")
    parser.add_argument('unittest_args', nargs='*')
    args = parser.parse_args()
    display = args.display
    verbose = args.verbose
    sys.argv[1:] = args.unittest_args
    lsst.utils.tests.init()
    unittest.main()
