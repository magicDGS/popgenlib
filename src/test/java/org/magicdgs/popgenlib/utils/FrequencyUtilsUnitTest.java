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

package org.magicdgs.popgenlib.utils;

import org.magicdgs.popgenlib.PopGenLibTest;

import org.apache.commons.math3.util.Pair;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class FrequencyUtilsUnitTest extends PopGenLibTest {

    @DataProvider(name = "goodFrequencies")
    public Object[][] getGoodFrequencies() {
        return new Object[][] {
                // only one value
                {Collections.singletonList(1d)},
                // with zero frequency
                {Arrays.asList(0d, 1d)},
                {Arrays.asList(1d, 0d)},
                {Arrays.asList(0.5, 0.5, 0d)},
                // two frequencies
                {Arrays.asList(0.5, 0.5)},
                {Arrays.asList(0.12, 0.88)},
                {Arrays.asList(0.225, 0.775)},
                {Arrays.asList(0.376789, 0.623211)},
                // three frequencies
                {Arrays.asList(1 / 3d, 1 / 3d, 1 / 3d)},
                {Arrays.asList(0.5, 0.25, 0.25)},
                {Arrays.asList(0.5, 0.255, 0.245)},
                // 10 frequencies
                {Arrays.asList(0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1)}
        };
    }

    @Test(dataProvider = "goodFrequencies")
    public void testValidateFrequencies(final List<Double> frequencies) {
        FrequencyUtils.validateFrequencies(frequencies);
    }

    @DataProvider(name = "badFrequencies")
    public Object[][] getBadFrequencies() {
        return new Object[][] {
                {null},
                {Collections.emptyList()},
                {Collections.singletonList(1.00001)},
                {Collections.singletonList(10d)},
                {Arrays.asList(1d, 1d)},
                {Arrays.asList(1 / 3d, 0.34)},
                {Collections.singletonList(null)},
                {Arrays.asList(0.5, null, 0.5)},
                {Arrays.asList(null, 0.5, 0.5)},
                {Arrays.asList(0.5, 0.5, null)},
                {Arrays.asList(1d, -1d)},
                {Arrays.asList(0d, 1d, 0.5, -0.5)}
        };
    }

    @Test(dataProvider = "badFrequencies", expectedExceptions = FrequencyUtils.IllegalFrequencyException.class)
    public void testInvalidValidateFrequencies(final List<Double> frequencies) {
        FrequencyUtils.validateFrequencies(frequencies);
    }

    @DataProvider
    public Object[][] frequenciesToSort() {
        return new Object[][] {
                // not modified
                {Collections.singletonList(1d), Collections.singletonList(1d)},
                {Arrays.asList(0.5, 0.5), Arrays.asList(0.5, 0.5)},
                // already sorted
                {Arrays.asList(0.7, 0.3), Arrays.asList(0.7, 0.3)},
                // sort them
                {Arrays.asList(0.3, 0.7), Arrays.asList(0.7, 0.3)}
        };
    }

    @Test(dataProvider = "frequenciesToSort")
    public void testSortFrequencies(final List<Double> frequencies, final List<Double> expected) {
        // make immutable to check that it is not modified inside
        final List<Double> immutable = Collections.unmodifiableList(frequencies);
        Assert.assertEquals(FrequencyUtils.sortFrequencies(immutable), expected);
    }

    @DataProvider(name = "badCounts")
    public Object[][] invalidCounts() {
        return new Object[][] {
                // null, empty or null-containing lists
                {null},
                {Collections.emptyList()},
                {Collections.singletonList(null)},
                // only zeroes
                {Collections.singletonList(0)},
                {Arrays.asList(0, 0)},
                // negative value
                {Arrays.asList(1, -1)}
        };
    }

    @Test(dataProvider = "badCounts", expectedExceptions = FrequencyUtils.IllegalFrequencyException.class)
    public void testInvalidCountsToFrequencies(final List<Integer> counts) {
        FrequencyUtils.countsToFrequencies(counts);
    }

    @DataProvider(name = "countsData")
    public Object[][] getCountsData() {
        return new Object[][] {
                // only one class
                {Collections.singletonList(1), 1, Collections.singletonList(1d)},
                {Collections.singletonList(10), 10, Collections.singletonList(1d)},
                // two classes
                {Arrays.asList(1, 1), 2, Arrays.asList(0.5, 0.5)},
                {Arrays.asList(1, 2), 3, Arrays.asList(1 / 3d, 2 / 3d)},
                {Arrays.asList(2, 1), 3, Arrays.asList(2 / 3d, 1 / 3d)},
                {Arrays.asList(1, 3), 4, Arrays.asList(0.25, 0.75)},
                // three classes
                {Arrays.asList(1, 1, 1), 3, Arrays.asList(1 / 3d, 1 / 3d, 1 / 3d)},
                {Arrays.asList(1, 2, 3), 6, Arrays.asList(1 / 6d, 2 / 6d, 3 / 6d)},
                // with zeroes
                {Arrays.asList(1, 0), 1, Arrays.asList(1d, 0d)},
                {Arrays.asList(10, 0), 10, Arrays.asList(1d, 0d)},
                {Arrays.asList(1, 1, 0), 2, Arrays.asList(0.5, 0.5, 0d)},
                {Arrays.asList(0, 20, 0, 90), 110, Arrays.asList(0.0, 20 / 110d, 0.0, 90 / 110d)}
        };
    }

    @Test(dataProvider = "countsData")
    public void testCountsToFrequencies(final List<Integer> counts, final int expectedCounts,
            final List<Double> expectedFrequencies) {
        // check if the result is as expected
        final Pair<Integer, List<Double>> result = FrequencyUtils.countsToFrequencies(counts);
        Assert.assertEquals(result.getFirst().intValue(), expectedCounts);

        final List<Double> freqs = result.getSecond();

        Assert.assertEquals(freqs.size(), expectedFrequencies.size(), "not equal length");
        for (int i = 0; i < expectedFrequencies.size(); i++) {
            // test with the maximum number of digits for double with precision
            Assert.assertEquals(freqs.get(i), expectedFrequencies.get(i), 1e-15,
                    "not equal value at " + i);
        }

        // and the validation of the frequencies should not blow up
        FrequencyUtils.validateFrequencies(result.getSecond());
    }
}