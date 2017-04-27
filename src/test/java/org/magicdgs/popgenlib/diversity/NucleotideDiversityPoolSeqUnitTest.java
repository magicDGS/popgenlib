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

import org.apache.commons.math3.distribution.BinomialDistribution;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class NucleotideDiversityPoolSeqUnitTest extends PopGenLibTest {


    @DataProvider(name = "piCounts")
    public Object[][] getAlleleCountsForTajimasPi() {
        // order of data
        // 1. minor allele count
        // 2. pool-size
        // 3. list of allele counts
        // 4. expected tajima's pi
        return new Object[][] {
                // all the following tests are ported from PoPoolation1:
                {2, 50, Arrays.asList(5, 5, 0, 0), 0.7147757},
                {2, 50, Arrays.asList(50, 50, 0, 0), 0.5187249},
                // pool-size 100, minor count 1
                {1, 100, Arrays.asList(5, 5, 0, 0), 0.5611672},
                {1, 100, Arrays.asList(6, 4, 0, 0), 0.5387205},
                {1, 100, Arrays.asList(70, 0, 30, 0), 0.4285277},
                {1, 100, Arrays.asList(0, 20, 0, 90), 0.303283},
                // Coverage influence
                // pool-size 100, minor count 2
                {2, 100, Arrays.asList(0, 0, 6, 4), 0.6858318},
                {2, 100, Arrays.asList(0, 0, 12, 8), 0.5648952},
                {2, 100, Arrays.asList(0, 0, 60, 40), 0.495659},
                // minor count influence
                {3, 100, Arrays.asList(0, 60, 0, 40), 0.5052888},
                {4, 100, Arrays.asList(0, 0, 60, 40), 0.5161168},
                // pool-size influence
                {2, 10, Arrays.asList(0, 60, 0, 40), 0.5387245},
                {2, 50, Arrays.asList(0, 60, 0, 40), 0.4979759},
                {1, 100, Arrays.asList(0, 20, 0, 480), 0.0777312}
        };
    }

    @Test(dataProvider = "piCounts")
    public void testTajimasPi(final int minCount, final int poolSize,
            final List<Integer> alleleCounts, final double expectedPi) throws Exception {
        Assert.assertEquals(NucleotideDiversityPoolSeq.tajimasPi(minCount, poolSize, alleleCounts),
                expectedPi, STATISTICAL_PRECISION);
    }

    @DataProvider(name = "thetaSingleSite")
    public static Object[][] getSingleSiteForTheta() {
        // order of data:
        // 1. minor allele count
        // 2. pool-size
        // 3. coverage
        // 4. expected watterson's theta
        return new Object[][] {
                // all the following tests are ported from PoPoolation1:
                // coverage influence
                {2, 100, 4, 2.0002},
                {2, 100, 10, 0.5822476},
                {2, 100, 20, 0.4010385},
                {2, 100, 100, 0.2423092},
                // minimum count influence
                {1, 100, 100, 0.2119774},
                {3, 100, 100, 0.2735117},
                {4, 100, 100, 0.3017763},
                // pool-size influence
                {2, 10, 100, 0.3535304},
                {2, 50, 100, 0.2489906},
                {2, 100, 500, 0.1946667}
        };
    }

    @DataProvider(name = "wrongPoolSeqParams")
    public Object[][] getIllegalArgsForSingleSiteTheta() {
        return new Object[][] {
                {0, 100},
                {1, 0}
        };
    }

    @Test(dataProvider = "thetaSingleSite")
    public void testSingleSiteWattersonsThetaCorrection(final int minCount, final int poolSize,
            final int coverage, final double expectedTheta) throws Exception {
        Assert.assertEquals(
                NucleotideDiversityPoolSeq.wattersonsTheta(1, minCount, poolSize, coverage),
                expectedTheta, STATISTICAL_PRECISION);

    }

    @Test(dataProvider = "thetaSingleSite")
    public void testTwoSitesWattersonsThetaCorrection(final int minCount, final int poolSize,
            final int coverage, final double expectedThetaForOneSite) {
        // test if it is the same as passing as a list of coverages
        Assert.assertEquals(
                NucleotideDiversityPoolSeq
                        .wattersonsTheta(minCount, poolSize, Collections.singletonList(coverage)),
                NucleotideDiversityPoolSeq.wattersonsTheta(1, minCount, poolSize, coverage));
        // test if it is the same as passing as a list of coverages with two sites
        Assert.assertEquals(
                NucleotideDiversityPoolSeq
                        .wattersonsTheta(minCount, poolSize, Arrays.asList(coverage, coverage)),
                NucleotideDiversityPoolSeq.wattersonsTheta(2, minCount, poolSize, coverage));
    }

    @DataProvider(name = "thetaSegregatingSites")
    public static Object[][] getSegregatingSitesForTheta() {
        // order of data:
        // 1. minor allele count
        // 2. pool-size
        // 3. list of coverages
        // 4. expected watterson's theta (computed with PoPoolation and multiply by the number of SNPs, because this is not an average)
        return new Object[][] {
                // tests computed with PoPoolation and a fake pileup
                // two segregating sites with same coverage
                {1, 10, Arrays.asList(20, 20), 0.372261772 * 2},
                // two segregating with different coverages
                {1, 10, Arrays.asList(20, 15), 0.380249799 * 2},
                {1, 10, Arrays.asList(15, 10), 0.406160223 * 2},
                // three segregating sites with different coverages
                {1, 10, Arrays.asList(20, 20, 15), 0.377587123 * 3},
                {1, 10, Arrays.asList(20, 15, 10), 0.394860739 * 3}
        };
    }


    @Test(dataProvider = "thetaSegregatingSites")
    public void testWattersonsThetaCorrection(final int minCount, final int poolSize,
            final List<Integer> coverages, final double expectedTheta) throws Exception {
        Assert.assertEquals(
                NucleotideDiversityPoolSeq.wattersonsTheta(minCount, poolSize, coverages),
                expectedTheta, STATISTICAL_PRECISION);
    }

    @DataProvider(name = "binomialParams")
    public Object[][] getBinomialParameters() {
        return new Object[][] {
                // higher coverage than pool-size
                {100, 500},
                // higher pool-size than coverage
                {500, 100},
                // equal pool-size and coverag
                {100, 100}
        };
    }

    @Test(dataProvider = "binomialParams")
    public void testBinomialProbability(final int poolSize, final int coverage) throws Exception {
        for (int i = 1; i < poolSize; i++) {
            for (int readCount = 0; readCount <= coverage; readCount++) {
                Assert.assertEquals(
                        NucleotideDiversityPoolSeq
                                .countProbability(readCount, coverage, poolSize, i),
                        new BinomialDistribution(coverage, (double) i / poolSize)
                                .probability(readCount),
                        STATISTICAL_PRECISION);
            }
        }
    }

    @Test(dataProvider = "wrongPoolSeqParams", expectedExceptions = IllegalArgumentException.class)
    public void testIllegalPoolSeqParams(final int minCount, final int poolSize) throws Exception {
        NucleotideDiversityPoolSeq.validatePoolSeqParams(minCount, poolSize);
    }

    @Test
    public void testIllegalParamsForWattersonsTheta() {
        // valid minCount and poolSize (this is tested in validatePoolSeqParams)
        final int minCount = 1;
        final int poolSize = 100;

        // empty list of coverages
        Assert.assertThrows(IllegalArgumentException.class,
                () -> NucleotideDiversityPoolSeq
                        .wattersonsTheta(minCount, poolSize, Collections.emptyList()));

        // wrong segregating sites
        Assert.assertThrows(IllegalArgumentException.class,
                () -> NucleotideDiversityPoolSeq.wattersonsTheta(-1, minCount, poolSize, 100));

        // wrong coverage
        Assert.assertThrows(IllegalArgumentException.class,
                () -> NucleotideDiversityPoolSeq.wattersonsTheta(1, minCount, poolSize, 0));

    }

}