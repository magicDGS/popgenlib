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

import org.magicdgs.popgenlib.diversity.NucleotideDiversity;
import org.magicdgs.popgenlib.utils.Verify;

import com.google.common.annotations.VisibleForTesting;
import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.util.Pair;

import java.util.function.Supplier;
import java.util.stream.IntStream;

/**
 * Static methods to compute Tajima's D.
 *
 * <p>REFERENCES:
 * <ul>
 *
 * <li><a href="http://www.genetics.org/content/123/3/585">
 * Tajima (1989): Statistical method for testing the neutral mutation hypothesis by DNA
 * polymorphism, Genetics 123(3).
 * </a></li>
 *
 * </ul>
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public final class TajimasD {

    // cannot be instantiated
    private TajimasD() {}

    /**
     * Computes Tajima's D using sample summaries using the formula 38 in
     * <a href="http://www.genetics.org/content/123/3/585">Tajima (1989)</a>
     *
     * @param numberOfSamples          number of samples used to get the number of segregating
     *                                 sites.
     * @param numberOfSegregatingSites number of segregating sites in the sample.
     * @param pairwiseDifferences      pair-wise differences between samples.
     */
    public static double tajimasD(final int numberOfSamples, final int numberOfSegregatingSites,
            final double pairwiseDifferences) {
        Verify.validate(numberOfSegregatingSites >= 0,
                () -> "Number of segregating sites should be a positive integer or 0: "
                        + numberOfSegregatingSites);
        Verify.validate(numberOfSamples > 1,
                () -> "numberOfSamples should be at least 2: " + numberOfSamples);

        // compute the variance constants
        final Pair<Double, Double> varianceConstants = varianceConstants(numberOfSamples);
        // compute the variance (denominator of formula 38)
        final double variance =
                FastMath.sqrt(varianceConstants.getFirst() * numberOfSegregatingSites
                        + varianceConstants.getSecond() * numberOfSegregatingSites * (
                        numberOfSegregatingSites - 1));
        // otherwise, follow formula 38
        return tajimasD(pairwiseDifferences,
                () -> NucleotideDiversity
                        .wattersonsTheta(numberOfSamples, numberOfSegregatingSites),
                variance);
    }

    // implementation of formula 38 with some advantages:
    // - checks for 0 variance
    // - computes theta only if the variance is not 0 (early optimization)
    @VisibleForTesting
    static double tajimasD(final double pairwiseDifferences, final Supplier<Double> thetaSupplier,
            final double variance) {
        // variance equals 0 returns directly
        if (variance == 0) {
            return 0;
        }
        // returns by computing theta
        return (pairwiseDifferences - thetaSupplier.get()) / variance;
    }

    /**
     * Computes the two constants for the variance of the difference between Tajima's &pi; and
     * Watterson's &theta; (d).
     *
     * <p>These values are used in the denominator of formula 38 from
     * <a href=http://www.genetics.org/content/123/3/585>Tajima (1989)</a> to compute the variance,
     * and corresponds to formulas 36 (e<SUB>1</SUB>) and 37 (e<SUB>2</SUB>).
     *
     * <p>These constants are based on Watterson's &theta; denominator, which is computed using and
     * approximated for large sample size.
     *
     * @param numberOfSamples number of samples to compute the denominator.
     *
     * @return e<SUB>1</SUB> and e<SUB>2</SUB>, in order.
     *
     * @see NucleotideDiversity#wattersonsDenominatorApproximation(int)
     */
    @VisibleForTesting
    static Pair<Double, Double> varianceConstants(final int numberOfSamples) {
        // Watterson's theta constants (formulas 3 and 4); a1 is computed with the approximation
        // a1 is watterson's denominator
        final double a1 = NucleotideDiversity.wattersonsDenominatorApproximation(numberOfSamples);
        final double a2 = IntStream.range(1, numberOfSamples)
                .mapToDouble(i -> 1d / (i * i)).sum();

        // Tajima's pi constants (formulas 8 and 9)
        final double b1 = (numberOfSamples + 1d)
                / (3d * (numberOfSamples - 1d));
        final double b2 = (2d * (FastMath.pow(numberOfSamples, 2) + numberOfSamples + 3d))
                / (9d * numberOfSamples * (numberOfSamples - 1d));

        // compute only once a1^2, because it is used twice
        final double a1Square = a1 * a1;

        // Variance constants (formulas 31 and 32)
        final double c1 = b1 - (1d / a1);
        final double c2 = b2
                - ((numberOfSamples + 2d) / (a1 * numberOfSamples))
                + (a2 / a1Square);

        // Variance estimation constants (formulas 36 and 37)
        final double e1 = c1 / a1;
        final double e2 = c2 / (a1Square + a2);

        // return as a pair
        return Pair.create(e1, e2);
    }
}
