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

import org.magicdgs.popgenlib.utils.FrequencyUtils;
import org.magicdgs.popgenlib.utils.Verify;

import com.google.common.annotations.VisibleForTesting;
import org.apache.commons.math3.special.Gamma;
import org.apache.commons.math3.util.FastMath;

import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Static methods to compute nucleotide diversity statistics.
 *
 * <p>REFERENCES:
 * <ul>
 *
 * <li><a href="http://www.genetics.org/content/123/3/585">
 * Tajima (1989): Statistical method for testing the neutral mutation hypothesis by DNA
 * polymorphism, Genetics 123(3).
 * </a></li>
 *
 * <li><a href="http://www.sciencedirect.com/science/article/pii/0040580975900209">
 * Watterson (1975): On the number of segregating sites in genetical models without recombination,
 * Theor. Popul. Biol. 7(2).
 * </a></li>
 *
 * </ul>
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public final class NucleotideDiversity {

    // cannot be instantiated
    private NucleotideDiversity() {}

    /**
     * Computes Tajima's &pi; for a single site using allele frequencies.
     *
     * <p>Corresponds to formula 12 in
     * <a href="http://www.genetics.org/content/123/3/585">Tajima (1989)</a>
     *
     * @param numberOfSamples   number of samples used for get the frequencies.
     * @param alleleFrequencies frequencies for each allele.
     */
    public static double tajimasPi(final int numberOfSamples,
            final List<Double> alleleFrequencies) {
        FrequencyUtils.validateFrequencies(alleleFrequencies);
        Verify.validate(numberOfSamples >= 2,
                () -> "numberOfSamples should be at least 2 for computing Tajima's Pi: "
                        + numberOfSamples);

        // computes the sum of p2
        final double pSquareSum =
                alleleFrequencies.stream().mapToDouble(i -> FastMath.pow(i, 2)).sum();
        // tajima's pi computed as in formula 12 in Tajima 1989
        return (numberOfSamples * (1 - pSquareSum)) / (numberOfSamples - 1);
    }

    /**
     * Computes Tajima's &pi; for a single site using allele counts.
     *
     * <p>Note: number of samples is computed using all allele counts.
     *
     * @param alleleCounts counts for each allele.
     *
     * @see #tajimasPi(int, List)
     */
    public static double tajimasPi(final List<Integer> alleleCounts) {
        Verify.nonEmpty(alleleCounts, () -> "allele counts");
        final int numberOfSamples = alleleCounts.stream().reduce(0, Integer::sum);
        final List<Double> freqs = alleleCounts.stream()
                .mapToDouble(c -> (double) c / (double) numberOfSamples)
                .boxed().collect(Collectors.toList());
        return tajimasPi(numberOfSamples, freqs);
    }

    /**
     * Computes Watterson's &theta; using the formula 1.4a in
     * <a href="http://www.sciencedirect.com/science/article/pii/0040580975900209">
     * Watterson (1975)</a>
     *
     * @param numberOfSamples          number of samples used to get the number of segregating
     *                                 sites.
     * @param numberOfSegregatingSites number of segregating sites in the sample.
     */
    public static double wattersonsTheta(final int numberOfSamples,
            final int numberOfSegregatingSites) {
        Verify.validate(numberOfSegregatingSites >= 0,
                () -> "Number of segregating sites should be a positive integer or 0: "
                        + numberOfSegregatingSites);
        Verify.validate(numberOfSamples > 1,
                () -> "numberOfSamples should be at least 2: " + numberOfSamples);
        // - S is the number of segregating sites
        // - a1 is the denominator of Watterson's theta (seee the method)
        // formula S / a1
        return numberOfSegregatingSites / wattersonsDenominatorApproximation(numberOfSamples);
    }

    /**
     * Gets the denominator of Watterson's &theta;
     * (formula 3 in <a href="http://www.genetics.org/content/123/3/585">Tajima (1989)</a>), which
     * is the {@code numberOfSamples}<SUP>th</SUP> -1 Harmonic Number.
     *
     * This method returns the exact value for up to 49 samples. Otherwise, it uses the
     * Euler–Mascheroni constant (&gamma;) and an approximation for the digamma distribution
     * (&psi;) for computing the n<SUP>th</SUP> Harmonic Number:
     *
     * <p>&gamma; + &psi;({@code numberOfSamples})
     *
     * @param numberOfSamples number of samples to compute the denominator.
     *
     * @throws IllegalArgumentException if the number of samples is lower than 2.
     * @see Gamma#digamma(double)
     */
    @VisibleForTesting
    static double wattersonsDenominatorApproximation(final int numberOfSamples) {
        // Defined as sum(1/j) for j=1 to j=n-1; where n is the number of samples
        // This is actually the (n-1)th Harmonic number (https://en.wikipedia.org/wiki/Harmonic_number)
        // This could be more efficiently computed using the formula implying the digamma distribution,
        // which is implemented in an approximated way in commons-math3 (good enough for our purposes).
        // In common-math3 they only use the fast algorithm if n is >= 49
        // thus, the limit for switching to the fast algorithm is 50 here,
        // because it will use the same number of iterations if not
        if (numberOfSamples < 49) {
            return IntStream.range(1, numberOfSamples).mapToDouble(i -> 1d / i).sum();
        }
        // The formula of the nth Harmonic number using the digamma function is defined as:
        // gamma + psi(n + 1); where gamma is the Euler-Mascheroni constant and psi the digamma function
        // because here we want the number numberOfSamples-1 Harmonic number, we can use directly
        // gamma + psi(numberOfSamples) if more than 50 (efficient computation)
        return Gamma.GAMMA + Gamma.digamma(numberOfSamples);
    }
}
