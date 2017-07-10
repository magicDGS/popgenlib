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
import org.apache.commons.math3.util.CombinatoricsUtils;
import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.util.Pair;

import java.util.List;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Static methods to compute nucleotide diversity statistics with corrections for Pool-Seq data.
 *
 * <p>REFERENCES:
 * <ul>
 *
 * <li><a href="http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0015925">
 * Kofler <i>et al.</i> (2010): PoPoolation: A Toolbox for Population Genetic Analysis of Next
 * Generation Sequencing Data from Pooled Individuals, PLOS ONE 6(1).
 * </a></li>
 *
 * </ul>
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public final class NucleotideDiversityPoolSeq {

    // cannot be instantiated
    private NucleotideDiversityPoolSeq() {}

    /**
     * Computes Tajima's &pi; for a single site using allele counts and correcting using the method
     * described in
     * <a href="http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0015925">
     * Kofler <i>et al.</i> (2010)</a> (page 7, first formula).
     *
     * @param minCount     minimum count used.
     * @param poolSize     number of samples pooled together.
     * @param alleleCounts counts for each allele.
     */
    public static double tajimasPi(final int minCount, final int poolSize,
            final List<Integer> alleleCounts) {
        validatePoolSeqParams(minCount, poolSize);
        // avoid recomputation of total coverage
        final Pair<Integer, List<Double>> freqs = FrequencyUtils.countsToFrequencies(alleleCounts);
        final double uncorrectedPi =
                NucleotideDiversity.tajimasPi(freqs.getFirst(), freqs.getSecond());
        final double correctionFactor = getPiCorrectionFactor(minCount, poolSize, freqs.getFirst());
        return correctionFactor * uncorrectedPi;
    }

    /**
     * Computes Watterson's &theta; using the corrected formula from
     * <a href="http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0015925">
     * Kofler <i>et al.</i> (2010)</a> (page 7, first formula).
     *
     * <p>Note: this implementation assumes the same minimum count, pool-size and coverage for all
     * the positions. For different coverages, use {@link #wattersonsTheta(int, int, List)}
     * instead.
     *
     * @param numberOfSegregatingSites number of segregating sites in the sample.
     * @param minCount                 minimum count used.
     * @param poolSize                 number of samples pooled together.
     * @param coverage                 non-zero coverage for segregating sites.
     *
     * @see #wattersonsTheta(int, int, List)
     */
    public static double wattersonsTheta(final int numberOfSegregatingSites, final int minCount,
            final int poolSize, final int coverage) {
        // verify params
        validatePoolSeqParams(minCount, poolSize);
        Verify.validate(numberOfSegregatingSites >= 0,
                () -> "Number of segregating sites should be a positive integer or 0: "
                        + numberOfSegregatingSites);
        Verify.validate(coverage > 0, () -> "Coverage should be a positive interger: " + coverage);


        // wattersons theta denominator is cancelled in the correction factor, so there is no need to compute it
        // this is the correction factor, using all posible read counts
        final double correctionFactor = IntStream.rangeClosed(minCount, coverage - minCount)
                .mapToDouble(readCount -> summationTerm(readCount, coverage, poolSize))
                .sum();
        // just divide by the correction factor, but not for the Watterson's denominator
        return numberOfSegregatingSites / correctionFactor;

    }

    /**
     * Computes Watterson's &theta; using the corrected formula from
     * <a href="http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0015925">
     * Kofler <i>et al.</i> (2010)</a> (page 7, second formula).
     *
     * <p>Note: the number of segregating sites is the size of the list with coverages.
     *
     * @param minCount                 minimum count used.
     * @param poolSize                 number of samples pooled together.
     * @param segregatingSitesCoverage non-zero coverages for segregating sites.
     *
     * @see #wattersonsTheta(int, int, int, int)
     */
    public static double wattersonsTheta(final int minCount, final int poolSize,
            final List<Integer> segregatingSitesCoverage) {
        Verify.nonEmpty(segregatingSitesCoverage, () -> "segregatingSitesCoverage");
        // for efficiency, all sites with the same coverage should compute the correction factor just once
        return segregatingSitesCoverage.stream()
                // categorize by the coverage
                .collect(Collectors.groupingBy(Function.identity(), Collectors.counting()))
                // for each entry
                .entrySet().stream()
                // the key is the coverage and the value the number of segregating sites in that category
                .mapToDouble(
                        entry -> wattersonsTheta(entry.getValue().intValue(), minCount, poolSize,
                                entry.getKey()))
                // sum them all afterwards
                .sum();
    }

    /**
     * Gets the correction factor for Tajima's &pi; defined in
     * <a href="http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0015925">
     * Kofler <i>et al.</i> (2010)</a> (page 7, first formula denominator).
     *
     * <p>The correction factor is based in the Pool-Size and coverage on the region, and setting a
     * minimum number of read counts to estimate frequencies.
     *
     * <p>Note: this should be applied to a Tajima's &pi; estimate already correcting by sample
     * size. This correction factor is <code>numberOfSamples - 1 / numberOfSamples</code>.
     *
     * @param minCount minimum count used to estimate frequencies.
     * @param poolSize number of individuals pooled together and used to estimate frequencies.
     * @param coverage coverage in the region used to estimate frequencies.
     *
     * @return correction factor, which should be multiplied to the uncorrected Tajima's &pi;.
     */
    private static double getPiCorrectionFactor(final int minCount, final int poolSize,
            final int coverage) {
        // for all posible read counts
        return 1d /
                IntStream.rangeClosed(minCount, coverage - minCount)
                        .mapToDouble(readCount -> {
                            // compute the first component, which is constant
                            final double term1 =
                                    2d * readCount * (coverage - readCount) / (coverage * (coverage
                                            - 1d));
                            // multiply this constant by the previous one
                            return term1 * summationTerm(readCount, coverage, poolSize);
                        })
                        .sum();
    }

    /**
     * Helper method to re-use in the correction factor for Tajima's &pi; and Watterson's &theta;.
     *
     * @return the second term in the divisor last formula in PoPoolation notes.
     */
    private static double summationTerm(final int readCount, final int coverage,
            final int poolSize) {
        // summation from k = 1 to k = poolSize - 1
        return IntStream.range(1, poolSize)
                // probability of each k over k
                // TODO: this could be substituted by:
                // TODO: probability = new BinomialDistribution(coverage, (double) k / poolSize).probability(readCount)
                // TODO: probability / k
                // TODO: needs tests for efficiency (see https://github.com/magicDGS/popgenlib/issues/23)
                .mapToDouble(k -> countProbability(readCount, coverage, poolSize, k) / k)
                .sum();
    }

    /**
     * Probability of having a first allele count of m in C reads from a pool of n with first
     * allele count of i (m is the allele count in the reads, i is teh allele count in the pool).
     * Corresponds to formula 2 in PoPoolation notes.
     *
     * <p>Note: this is the probability of {@code readCount} assuming a binomial distribution with
     * parameters n={@code coverage} and p={@code poolCount / poolSize}.
     *
     * @return binomial probability.
     */
    @VisibleForTesting
    static double countProbability(final int readCount, final int coverage, final int poolSize,
            final int poolCount) {
        final double coef = CombinatoricsUtils.binomialCoefficientDouble(coverage, readCount);
        final double firstTerm = FastMath.pow((double) poolCount / poolSize, readCount);
        final double seconfTerm =
                FastMath.pow((double) (poolSize - poolCount) / poolSize, (coverage - readCount));

        return coef * firstTerm * seconfTerm;
    }

    /**
     * Validates the Pool-Seq parameters for the correction factor (minimum read count and
     * pool-size).
     */
    @VisibleForTesting
    static final void validatePoolSeqParams(final int minCount, final int poolSize) {
        Verify.validate(minCount > 0,
                () -> "Minimum count should be a positive integer: " + minCount);
        Verify.validate(poolSize > 0, () -> "Pool-Size should be a positive integer: " + poolSize);
    }

}
