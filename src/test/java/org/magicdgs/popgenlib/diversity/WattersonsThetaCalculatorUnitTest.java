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

package org.magicdgs.popgenlib.diversity;

import org.magicdgs.popgenlib.PopGenLibTest;

import com.google.common.cache.CacheBuilder;
import com.google.common.cache.LoadingCache;
import org.testng.Assert;
import org.testng.annotations.Test;

/**
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class WattersonsThetaCalculatorUnitTest extends PopGenLibTest {

    // to test the calculation
    private final static WattersonsThetaCalculator calculator = new WattersonsThetaCalculator();

    @Test(dataProvider = "thetaSegregatingSites", dataProviderClass = NucleotideDiversityUnitTest.class)
    public void testWattersonsTheta(final int numberOfSamples, final int numberOfSegregatingSites,
            final double expectedTheta) throws Exception {
        // test to compute it twice and still it returns the same value
        Assert.assertEquals(calculator.wattersonsTheta(numberOfSamples, numberOfSegregatingSites),
                expectedTheta, STATISTICAL_PRECISION);
        Assert.assertEquals(calculator.wattersonsTheta(numberOfSamples, numberOfSegregatingSites),
                expectedTheta, STATISTICAL_PRECISION);
    }

    @Test(dataProvider = "thetaInvalidParams", dataProviderClass = NucleotideDiversityUnitTest.class, expectedExceptions = IllegalArgumentException.class)
    public void testWattersonsThetaInvalidValues(final int numberOfSamples,
            final int segregatingSite) {
        calculator.wattersonsTheta(numberOfSamples, segregatingSite);
    }

    @Test
    public void testValuesAreCached() throws Exception {
        // init params
        final WattersonsThetaCalculator calculator =
                new WattersonsThetaCalculator(CacheBuilder.newBuilder().recordStats());
        int hits = 0;
        int missing = 0;
        // nothing at the beginning
        testCacheStats(calculator.denominatorCache, missing, hits);

        // compute theta 10 times with the same params (new missing)
        missing++;
        for (; hits <= 10; hits++) {
            calculator.wattersonsTheta(10, 2);
            testCacheStats(calculator.denominatorCache, missing, hits);
        }
        // last hit should be removed because it is not true
        hits--;

        // include a new value in the cache (new missing)
        calculator.wattersonsTheta(100, 2);
        missing++;
        testCacheStats(calculator.denominatorCache, missing, hits);
        // same number of samples (new hit, different segregating sites should not affect)
        calculator.wattersonsTheta(100, 10);
        hits++;
        testCacheStats(calculator.denominatorCache, missing, hits);

    }

    private static final void testCacheStats(final LoadingCache cache, final int missing,
            final int hits) throws Exception {
        Assert.assertEquals(cache.stats().hitCount(), hits);
        Assert.assertEquals(cache.stats().missCount(), missing);
        // size should be the same as missing for this test
        Assert.assertEquals(cache.size(), missing);
    }

}