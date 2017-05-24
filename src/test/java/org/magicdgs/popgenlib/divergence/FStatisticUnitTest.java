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

package org.magicdgs.popgenlib.divergence;

import org.magicdgs.popgenlib.PopGenLibTest;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.List;

/**
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class FStatisticUnitTest extends PopGenLibTest {

    @DataProvider(name = "pairwiseData")
    public Object[][] getPairwiseData() {
        return new Object[][] {
                // edge case 1: equal frequencies
                {10, Arrays.asList(0.1, 0.9), 10, Arrays.asList(0.1, 0.9), 0d},
                // edge case 2: fixed in both populations
                {10, Arrays.asList(1d, 0d), 10, Arrays.asList(0d, 1d), 1d},
                // small differences TODO: failing!
                // {30, Arrays.asList(0.33333333, 0.66666667), 50, Arrays.asList(0.4, 0.6), 0.01181525}
        };
    }

    @Test(dataProvider = "pairwiseData")
    public void testPairwiseFst(final int numberOfSamples1, final List<Double> alleleFrequencies1,
            final int numberOfSamples2, final List<Double> alleleFrequencies2,
            final double expectedFst)
            throws Exception {
        Assert.assertEquals(FStatistic.pairwiseFst(
                numberOfSamples1, alleleFrequencies1, numberOfSamples2, alleleFrequencies2),
                expectedFst, STATISTICAL_PRECISION);
    }

}