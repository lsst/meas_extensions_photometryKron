from builtins import str
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
import unittest
import uuid

import lsst.utils.tests
import lsst.afw.table as afwTable
import lsst.meas.extensions.photometryKron as photKron
from lsst.daf.base import PropertyList


def getTableDelimeter(schema):
    """
    Return the character(s) used as delimeters within the given schema.

    Should be "_" in the current afw.table implementation. This functionality
    is not otherwise exposed to Python.
    """
    a, b = str(uuid.uuid1()), str(uuid.uuid1())
    return schema.join(a, b).replace(a, '').replace(b, '')


class KronFlagHandlerTestCase(unittest.TestCase):
    """Test the FlagHandler used for Kron photometry"""

    def testFlagDefinitions(self):
        """
        Check flag order.

        Flags must be added to the flag handler in the same order that they
        are defined in the algorithm.
        """
        control = photKron.KronFluxControl()
        name = "kronTest"
        schema = afwTable.SourceTable.makeMinimalSchema()
        algMeta = PropertyList()

        # Add the output fields -- including flags -- to the schema.
        photKron.KronFluxAlgorithm(control, name, schema, algMeta)

        # Fetch a list of all flag fields, in the order that they were added
        # to the schema (and hence the order they were added to the FlagHandler)
        flagFieldNames = schema.extract("%s_flag*" % (name,), ordered=True).keys()
        # Iterate over each flag field, checking that they were enumerated in
        # the algorithm in the same order as in the FlagHandler.
        for i, flagFieldName in enumerate(flagFieldNames):
            if flagFieldName == "%s_flag" % (name,):
                # The generic "failure" flag is written into the schema as $name_flag.
                self.assertEqual(i, photKron.KronFluxAlgorithm.FAILURE.number)
            else:
                # Other flags are referenced by name. We assert that the
                # enumeration name (e.g. BAD_RADIUS) is an upper-case version
                # of the schema field name (e.g. flag_bad_radius), with the
                # "flag_" prefix stripped.
                flagName = flagFieldName.split(getTableDelimeter(schema), 2)[-1]
                self.assertEqual(i, getattr(photKron.KronFluxAlgorithm, flagName.upper()).number)

        # Check that the number of enumerated flags matches the number of flag
        # fields in the schema.
        self.assertEqual(len(photKron.KronFluxAlgorithm.getFlagDefinitions()), len(flagFieldNames))


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
