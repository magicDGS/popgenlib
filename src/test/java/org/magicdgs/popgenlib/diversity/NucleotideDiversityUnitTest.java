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

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.stream.IntStream;

/**
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class NucleotideDiversityUnitTest extends PopGenLibTest {

    @DataProvider(name = "piFrequencies")
    public Object[][] getAlleleFrequenciesForTajimasPi() {
        // order of data:
        // 1. number of samples
        // 2. allele frequencies
        // 3. Expected Tajima's pi
        return new Object[][] {
                // low diversity
                {100, Arrays.asList(0.95, 0.05), 0.0959596},
                {50, Arrays.asList(0.95, 0.05), 0.09693878},
                {10, Arrays.asList(0.95, 0.05), 0.1055556},
                // intermediate frequency
                {100, Arrays.asList(0.1, 0.9), 0.1818182},
                {50, Arrays.asList(0.1, 0.9), 0.1836735},
                {10, Arrays.asList(0.1, 0.9), 0.2000000},
                // high diversity
                {100, Arrays.asList(0.55, 0.45), 0.5000000},
                {50, Arrays.asList(0.55, 0.45), 0.505102},
                {10, Arrays.asList(0.55, 0.45), 0.55000000},
                // frequency ordered in the reverse way produce the same result
                {100, Arrays.asList(0.45, 0.55), 0.5000000},
                {50, Arrays.asList(0.45, 0.55), 0.505102},
                {10, Arrays.asList(0.45, 0.55), 0.55000000},
                // both allele frequencies are equal
                {100, Arrays.asList(0.5, 0.5), 0.5050505},
                {50, Arrays.asList(0.5, 0.5), 0.5102041},
                {10, Arrays.asList(0.5, 0.5), 0.5555555},
                // monomorphic site has pi = 0 independently of the number of samples
                // or if alleles with frequency zero are passed
                {10, Collections.singletonList(1d), 0d},
                {50, Collections.singletonList(1d), 0d},
                {100, Arrays.asList(1d, 0d, 0d), 0d},
                // tri-allelic site
                {10, Arrays.asList(0.5, 0.4, 0.1), 0.6444444},
                {50, Arrays.asList(0.5, 0.4, 0.1), 0.5918367},
                {100, Arrays.asList(0.5, 0.4, 0.1), 0.5858586}
        };
    }

    @Test(dataProvider = "piFrequencies")
    public void testTajimasPi(final int numberOfSamples, final List<Double> alleleFrequencies,
            final double expectedPi) {
        final double result = NucleotideDiversity.tajimasPi(numberOfSamples, alleleFrequencies);
        Assert.assertEquals(result, expectedPi, STATISTICAL_PRECISION, "freqs");
    }

    @DataProvider(name = "piCounts")
    public Object[][] getAlleleCountsForTajimasPi() {
        // this is less exhaustive than the one with allele frequencies
        // this is because it is exactly the same method, but computing the allele frequencies
        return new Object[][] {
                // monomorphic site with two samples
                {Collections.singletonList(2), 0d},
                {Arrays.asList(2, 0, 0), 0d},
                // same counts
                {Arrays.asList(5, 5), 0.5555555},
                // tri-allelic
                {Arrays.asList(5, 4, 1), 0.6444444}
        };
    }

    @Test(dataProvider = "piCounts")
    public void testTajimasPiWithCounts(final List<Integer> alleleCounts, final double expextedPi) {
        Assert.assertEquals(NucleotideDiversity.tajimasPi(alleleCounts), expextedPi,
                STATISTICAL_PRECISION);
    }

    @DataProvider(name = "piInvalidParamsWithFrequencies")
    public Object[][] getInvalidParamsForTajimasPiWithFrequencies() {
        final List<Double> alleleFreqs = Collections.singletonList(1d);
        return new Object[][] {
                // testing only invalid number of samples
                {-1, alleleFreqs},
                {0, alleleFreqs},
                {1, alleleFreqs}
        };
    }

    @Test(dataProvider = "piInvalidParamsWithFrequencies", expectedExceptions = IllegalArgumentException.class)
    public void testInvalidTajimasPiWithFrequencies(final int numberOfSamples,
            final List<Double> alleleFrequencies) throws Exception {
        NucleotideDiversity.tajimasPi(numberOfSamples, alleleFrequencies);
    }


    @DataProvider(name = "thetaSegregatingSites")
    public static Object[][] getSegregatingSitesForWattersonsTheta() {
        // order of data:
        // 1. number of samples
        // 2. number of segregating sites
        // 3. Expected Watterson's Theta
        return new Object[][] {
                // only 2 samples and 0 segregating sites (edge case)
                {2, 0, 0d},
                // only 2 samples and 1 segregating site (edge case)
                {2, 1, 1d},
                // for non-segregating sites is always 0 independently of the samples
                {100, 0, 0d},
                {10, 0, 0d},
                // for only one segregating and different number of samples
                {20, 1, 0.2818696},
                {100, 1, 0.1931480},
                // two segregating sites
                {20, 2, 0.5637392},
                {100, 2, 0.3862960},
                // three segregating sites
                {20, 3, 0.8456088},
                {100, 3, 0.5794439},
                // same number of segregating sites and samples
                {20, 20, 5.6373922},
                {100, 100, 19.3147978},
                // very large sample size to test the approximation
                {10000, 3, 0.3065132},
                {10000, 100, 10.2171074},
                {10000, 200, 20.4342147},
                {10000, 1000, 102.1710736}
        };
    }

    @Test(dataProvider = "thetaSegregatingSites")
    public void testWattersonsTheta(final int numberOfSamples, final int numberOfSegregatingSites,
            final double expectedTheta) throws Exception {
        Assert.assertEquals(
                NucleotideDiversity.wattersonsTheta(numberOfSamples, numberOfSegregatingSites),
                expectedTheta, STATISTICAL_PRECISION);
    }

    @DataProvider(name = "thetaInvalidParams")
    public static Object[][] getInvalidParamsForWattersonsTheta() {
        return new Object[][] {
                // invalid number of samples
                {-1, 10}, {0, 10}, {1, 10},
                // invalid segregating sites
                {10, -1}
        };
    }

    @Test(dataProvider = "thetaInvalidParams", expectedExceptions = IllegalArgumentException.class)
    public void testWattersonsThetaInvalidValues(final int numberOfSamples,
            final int segregatingSites) throws Exception {
        NucleotideDiversity.wattersonsTheta(numberOfSamples, segregatingSites);
    }

    @DataProvider(name = "denominatorData")
    public Iterator<Object[]> samplesForDenominator() {
        return IntStream.concat(
                // small number of samples (1 will never be passed in the denominator)
                IntStream.range(2, 100),
                // large number of samples (important for Pool-Seq, for instance)
                // and it is where we gain more in performance
                IntStream.range(1000, 1100)
        ).mapToObj(i -> new Object[] {i}).iterator();
    }

    @Test(dataProvider = "denominatorData")
    public void testWattersonsThetaDenominator(final int numberOfSamples) {
        // computes as a summation to be sure that it is correct
        // this is very inefficient for large samples, that's why we use the approximation
        // in the main code, but the test should include it
        final double summation = IntStream.range(1, numberOfSamples).mapToDouble(i -> 1d / i).sum();
        Assert.assertEquals(NucleotideDiversity.wattersonsDenominatorApproximation(numberOfSamples),
                summation, STATISTICAL_PRECISION);

    }
}