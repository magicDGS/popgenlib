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

package org.magicdgs.popgenlib.neutralitytest;

import org.magicdgs.popgenlib.PopGenLibTest;

import org.apache.commons.math3.util.Pair;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

/**
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class TajimasDUnitTest extends PopGenLibTest {

    @DataProvider(name = "parametersForTajimasD")
    public Object[][] getParametersForTajimasD() {
        // order of data:
        // 1. number of samples
        // 2. number of segregating sites
        // 3. number of pair-wise differences
        // 4. Expected Tajima's D
        return new Object[][] {
                // examples from https://ocw.mit.edu/courses/health-sciences-and-technology/hst-508-quantitative-genomics-fall-2005/study-materials/tajimad1.pdf
                {10, 16, 3.888889, -1.446172},
                {77, 103, 8.438483, -2.021749},
                {72, 88, 15.339984, -0.5258013}
        };
    }

    @Test(dataProvider = "parametersForTajimasD")
    public void testTajimasD(final int numberOfSamples, final int numberOfSegregatingSites,
            final double pairwiseDifferences,
            final double expectedTajimasD) {
        Assert.assertEquals(
                TajimasD.tajimasD(numberOfSamples, numberOfSegregatingSites, pairwiseDifferences),
                expectedTajimasD, STATISTICAL_PRECISION);
    }

    @Test
    public void testTajimasDWithZeroVariance() {
        // simple fake test because it is difficult to get an example with variance 0
        // but we have to protect to this case
        Assert.assertEquals(TajimasD.tajimasD(0, () -> 0d, 0), 0d);
    }

    @DataProvider(name = "varianceConstants")
    public Object[][] getVarianceConstantsData() {
        // order of data:
        // 1. number of samples
        // 2. expected e1
        // 3. expected e2
        return new Object[][] {
                // examples from https://ocw.mit.edu/courses/health-sciences-and-technology/hst-508-quantitative-genomics-fall-2005/study-materials/tajimad1.pdf
                {10, 0.0190605, 0.004949},
                {77, 0.0282075, 0.0033735},
                {72, 0.028143, 0.0034226}
        };
    }

    @Test(dataProvider = "varianceConstants")
    public void testVarianceConstants(final int numberOfSamples,
            final double expectedE1, final double expectedE2) {
        final Pair<Double, Double> constants = TajimasD.varianceConstants(numberOfSamples);
        Assert.assertEquals(constants.getFirst(), expectedE1, STATISTICAL_PRECISION, "e1");
        Assert.assertEquals(constants.getSecond(), expectedE2, STATISTICAL_PRECISION, "e2");
    }


    @DataProvider(name = "invalidParamsForTajimasD")
    public Object[][] getInvalidParamsForTajimasD() {
        return new Object[][] {
                // invalid number of samples
                {-1, 100, 1.11},
                {0, 100, 1.11},
                {1, 100, 1.11},
                // invalid number of segregating sites
                {2, -1, 1.11},
                {2, -2, 1.11}
        };
    }

    @Test(dataProvider = "invalidParamsForTajimasD", expectedExceptions = IllegalArgumentException.class)
    public void testInvalidTajimasD(final int numberOfSamples, final int numberOfSegregatingSites,
            final double pairwiseDifferences) throws Exception {
        TajimasD.tajimasD(numberOfSamples, numberOfSegregatingSites, pairwiseDifferences);
    }

}