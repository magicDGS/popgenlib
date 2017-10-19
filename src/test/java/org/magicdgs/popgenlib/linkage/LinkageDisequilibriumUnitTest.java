/*
 * MIT License
 *
 * Copyright (c) 2017 Daniel Gomez-Sanchez
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

package org.magicdgs.popgenlib.linkage;

import org.magicdgs.popgenlib.PopGenLibTest;
import org.magicdgs.popgenlib.utils.FrequencyUtils;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

/**
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class LinkageDisequilibriumUnitTest extends PopGenLibTest {

    @DataProvider(name = "dData")
    public Object[][] getDdata() {
        return new Object[][] {
                // equilibrium
                {0.5, 0.5, 0.25, 0},
                {0.7, 0.8, 0.7 * 0.8, 0},
                // independent of the input order, the same D (positive
                {0.7, 0.8, 0.6, 0.04},
                {0.8, 0.7, 0.6, 0.04},
                // negative D
                {0.5, 0.5, 0.1, -0.15},
                // complete disequilibrium
                {0.5, 0.5, 0, -0.25},
                {0.5, 0.5, 0.5, 0.25}
        };
    }

    @Test(dataProvider = "dData")
    public void tesD(final double pA, final double pB, final double pAB, final double expectedD) {
        Assert.assertEquals(LinkageDisequilibrium.d(pA, pB, pAB), expectedD, STATISTICAL_PRECISION);
    }

    @DataProvider(name = "invalidDdata")
    public Object[][] getInvalidDdata() {
        return new Object[][] {
                {0.4, 0.9, 0.9},
                {0.9, 0.4, 0.9},
                {0.9, 0.9, 10},
                {0.9, 0.9, -1}
        };
    }

    @Test(dataProvider = "invalidDdata", expectedExceptions = FrequencyUtils.IllegalFrequencyException.class)
    public void testInvalidDargs(final double pA, final double pB, final double pAB) {
        LinkageDisequilibrium.d(pA, pB, pAB);
    }

    @DataProvider(name = "dPrimeData")
    public Object[][] getDprimeData() {
        return new Object[][] {
                // equilibrium
                {0.5, 0.5, 0.25, 0},
                {0.7, 0.8, 0.7 * 0.8, 0},
                // independent of the input order, the same D' (positive)
                {0.7, 0.8, 0.6, 0.2857143},
                {0.8, 0.7, 0.6, 0.2857143},
                // negative D'
                {0.5, 0.5, 0.1, -0.6},
                // complete disequilibrium
                {0.5, 0.5, 0, -1},
                {0.5, 0.5, 0.5, 1}
        };
    }

    @Test(dataProvider = "dPrimeData")
    public void tesDPrime(final double pA, final double pB, final double pAB,
            final double expectedD) {
        Assert.assertEquals(LinkageDisequilibrium.dPrime(pA, pB, pAB), expectedD,
                STATISTICAL_PRECISION);
    }

    @DataProvider(name = "significantR2")
    public Object[][] significatR2data() {
        return new Object[][] {
                {0.25, 36, 0.95, true},
                {0.10, 36, 0.95, false}
        };
    }

    @Test(dataProvider = "significantR2")
    public void testR2significantTest(final double r2, final int numberOfSamples,
            final double quantile, final boolean significant) {
        Assert.assertEquals(LinkageDisequilibrium.r2SignificantTest(r2, numberOfSamples, quantile),
                significant);
    }

    @DataProvider(name = "significantTest")
    public Object[][] getSignificantThresholdData() {
        // computed with R
        return new Object[][] {
                {28, 0.95, 0.137195},
                {36, 0.95, 0.1067072},
                {205, 0.95, 0.01873882},
                {28, 0.99, 0.2369606},
                {36, 0.99, 0.1843027},
                {205, 0.99, 0.03236535}
        };
    }

    @Test(dataProvider = "significantTest")
    public void testR2SignificantThreshold(final int nHaplotypes, final double quantile,
            final double expected) throws Exception {
        Assert.assertEquals(LinkageDisequilibrium.r2SignificantThreshold(nHaplotypes, quantile),
                expected, STATISTICAL_PRECISION);
    }

    @DataProvider
    public Object[][] invalidSignificantR2() {
        return new Object[][] {
                // invalid r2
                {-1, 100},
                {20, 100},
                // invalid numberOfSamples
                {0.5, 0},
                {0.5, 1}
        };
    }

    @Test(dataProvider = "invalidSignificantR2", expectedExceptions = IllegalArgumentException.class)
    public void invalidParamsForR2SignificantTest(final double r2, final int numberOfSamples) {
        LinkageDisequilibrium.r2SignificantTest(r2, numberOfSamples, 0.95);
    }

    @DataProvider(name = "rwData")
    public Object[][] getRwData() {
        return new Object[][] {
                // completely linked
                {7d / 11, 7d / 11, 7d / 11, 1d},
                // positive rw
                {5d / 11, 7d / 11, 7d / 11, 0.2142857},
                {5d / 11, 7d / 11, 6d / 11, 0.4485426},
                // negative rw
                {3d / 11, 7d / 11, 6d / 11, -0.3105295},
                // with singletons
                {7d / 11, 7d / 11, 10d / 11, 0.41833},
                {6d / 11, 7d / 11, 10d / 11, -0.2390457},
                {5d / 11, 6d / 11, 10d / 11, -0.2886751}

        };
    }

    @Test(dataProvider = "rwData")
    public void testRw(final double pAB, final double pA, final double pB, final double expected)
            throws Exception {
        Assert.assertEquals(LinkageDisequilibrium.rw(pA, pB, pAB), expected, STATISTICAL_PRECISION);
    }

    @DataProvider(name = "r2Data")
    public Object[][] getR2Data() {
        return new Object[][] {
                // completely linked
                {7d / 11, 7d / 11, 7d / 11, 1d},
                // positive rw
                {5d / 11, 7d / 11, 7d / 11, 0.04591837},
                {5d / 11, 7d / 11, 6d / 11, 0.2011905},
                // negative rw
                {3d / 11, 7d / 11, 6d / 11, 0.09642857},
                // with singletons
                {7d / 11, 7d / 11, 10d / 11, 0.175},
                {6d / 11, 7d / 11, 10d / 11, 0.05714286},
                {5d / 11, 6d / 11, 10d / 11, 0.08333333}

        };
    }


    @Test(dataProvider = "r2Data")
    public void testR2(final double pAB, final double pA, final double pB, final double expected)
            throws Exception {
        Assert.assertEquals(LinkageDisequilibrium.r2(pA, pB, pAB), expected, STATISTICAL_PRECISION);
    }

    @DataProvider(name = "maxR2data")
    public Object[][] getMaxR2Data() {
        return new Object[][] {
                // equal frequencies
                {7d / 11, 7d / 11, 1d},
                // sligthly different frequencies
                {7d / 11, 6d / 11, 0.6857143},
                {6d / 11, 7d / 11, 0.6857143},
                // with singletons
                {7d / 11, 10d / 11, 0.175},
                {6d / 11, 10d / 11, 0.12}

        };
    }

    @Test(dataProvider = "maxR2data")
    public void testMaxR2(final double pA, final double pB, final double expected)
            throws Exception {
        Assert.assertEquals(LinkageDisequilibrium.r2max(pA, pB), expected, STATISTICAL_PRECISION);
    }

    @DataProvider(name = "invalidMaxR2")
    public Object[][] getInvalidMaxR2args() {
        return new Object[][] {
                {-1, 1},
                {0, 1},
                {0.2, 1},
                {1, -1},
                {1, 0},
                {1, 0.2},
        };
    }

    @Test(dataProvider = "invalidMaxR2", expectedExceptions = IllegalArgumentException.class)
    public void invalidMaxR2args(final double pA, final double pB) {
        LinkageDisequilibrium.r2max(pA, pB);
    }

    @DataProvider(name = "r2PrimeData")
    public Object[][] getR2PrimeData() {
        return new Object[][] {
                // completely linked
                {7d / 11, 7d / 11, 7d / 11, 1d},
                // positive rw
                {5d / 11, 7d / 11, 7d / 11, 0.04591837},
                {5d / 11, 7d / 11, 6d / 11, 0.2011905 / 0.6857143},
                // negative rw
                {3d / 11, 7d / 11, 6d / 11, 0.09642857 / 0.6857143},
                // with singletons
                {7d / 11, 7d / 11, 10d / 11, 0.175 / 0.175},
                {6d / 11, 7d / 11, 10d / 11, 0.05714286 / 0.175},
                {5d / 11, 6d / 11, 10d / 11, 0.08333333 / 0.12}

        };
    }

    @Test(dataProvider = "r2PrimeData")
    public void testR2Prime(final double pAB, final double pA, final double pB,
            final double expected) throws Exception {
        Assert.assertEquals(LinkageDisequilibrium.r2Prime(pA, pB, pAB), expected,
                STATISTICAL_PRECISION);
    }

}