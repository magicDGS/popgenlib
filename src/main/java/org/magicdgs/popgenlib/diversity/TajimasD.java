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

import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.util.Pair;

import java.util.function.Supplier;
import java.util.stream.IntStream;

/**
 * TODO: tajimas D documentation
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class TajimasD {



    // TODO: document
    public static double tajimasD(final int numberOfSamples, final int piEstimate, final int segregatingSites) {
        return tajimasD(piEstimate, segregatingSites,
                () -> NucleotideDiversity.wattersonsTheta(numberOfSamples, segregatingSites),
                varianceConstants(numberOfSamples));
    }

    // TODO: document
    // this is for allowing caching
    static double tajimasD(final int pi, final int S,
            Supplier<Double> wattersonsThetaSupplier, final Pair<Double, Double> varianceConstants) {
        // denominator formula 38 (Tajima 1989)
        final double variance = FastMath.sqrt(
                varianceConstants.getFirst() * S
                        + varianceConstants.getSecond() * S * (S - 1));

        // avoid division by 0; if not, formula 38
        return (variance == 0)
                ? 0
                : (pi - wattersonsThetaSupplier.get()) / variance;
    }

    /**
     * Computes the constants (e<SUB>1</SUB> and e<SUB>2</SUB>) for the variance of d (the
     * difference between Tajima's &pi; and Watterson's &theta;).
     *
     * These values are used in the denominator of formula 38 from
     * <a href=http://www.genetics.org/content/123/3/585>Tajima (1989)</a> to compute the variance.
     *
     * <p>TODO: write formula
     *
     * @param numberOfSamples number of samples to compute the denominator.
     */
    // TODO: this value could be cached
    static Pair<Double, Double> varianceConstants(final int numberOfSamples) {
        // all this values are defined in Tajima (1989)

        // WATERSON'S THETA CONSTANT (formulas 3)
        // a1 is watterson's denominator
        final double a1 = NucleotideDiversity.wattersonsDenominator(numberOfSamples);

        return varianceConstants(a1, numberOfSamples);
    }


    // this is necessary for allow caching of a1
    // TODO
    static Pair<Double, Double> varianceConstants(final double a1, final int numberOfSamples) {
        // all this values are defined in Tajima (1989)

        // WATTERSON'S THETA CONSTANT (formula 4)
        final double a2 = IntStream.range(1, numberOfSamples - 1)
                .mapToDouble(i -> 1d / FastMath.pow(i, 2)).sum();

        // TAJIMA'S PI CONSTANTS (formulas 8 and 9)
        final double b1 = (numberOfSamples + 1d)
                / (3d * (numberOfSamples - 1d));
        final double b2 = (2d * (FastMath.pow(numberOfSamples, 2) + numberOfSamples + 3d))
                / (9d * numberOfSamples * (numberOfSamples - 1d));

        // a1^2 is used twice, and square a value is demanding
        final double a1Square = FastMath.pow(a1, 2);

        // VARIANCE CONSTANTS (formulas 31 and 32)
        final double c1 = b1 - (1d / a1);
        final double c2 = b2
                - ((numberOfSamples + 2d) / (a1 * numberOfSamples))
                + (a2 / a1Square);

        // VARIANCE ESTIMATION CONSTANTS (formulas 36 and 37)
        final double e1 = c1 / a1;
        final double e2 = c2 / (a1Square + a2);

        return Pair.create(e1, e2);
    }
}
