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

package org.magicdgs.popgenlib.linkage;

import org.magicdgs.popgenlib.utils.FrequencyUtils;
import org.magicdgs.popgenlib.utils.Verify;

import org.apache.commons.math3.distribution.ChiSquaredDistribution;
import org.apache.commons.math3.util.FastMath;

/**
 * Static methods to compute linkage disequilibrium statistics.
 *
 * <p>Caveats:
 *
 * <ul>
 * <li>Methods in this class assume haplotype frequencies and requires phased data</li>
 * <li>Allele frequencies are assumed to be polarize, and only major frequencies provided.
 * This simplifies the formulas, and requires users to provide {@code 1 - p} if
 * {@code p > 0.5}.</li>
 * </ul>
 *
 *
 * <p>REFERENCES:
 *
 * <ul>
 *
 * <li><a href="http://www.genetics.org/content/49/1/49">
 * Lewontin (1964): The interaction of selection and linkage. I. General considerations; heterotic
 * models, Genetics 49.
 * </a></li>
 *
 * <li><a href="https://link.springer.com/article/10.1007%2FBF01245622">
 * Hill &amp; Robertson (1968): Linkage disequilibrium in finite populations, Theoretical and
 * Applied Genetics 38.
 * </a></li>
 *
 * <li><a href="http://www.genetics.org/content/78/3/937.long">
 * Langley &amp; Crow (1974): The direction of linkage disequilibrium, Genetics 78.
 * </a></li>
 *
 * <li><a href="http://www.sciencedirect.com/science/article/pii/S0040580908000609">
 * VanLiere &amp; Rosenberg (2008): Mathematical properties of the <i>r</i><SUP>2</SUP> measure of
 * linkage disequilibrium, Theoretical Population Biology 74(1).
 * </a></li>
 *
 * <li><a href="https://books.google.com/books/about/Elements_of_Evolutionary_Genetics.html?id=dgNFAQAAIAAJ">
 * Charlesworth &amp; Charlesworth (2012): Multiple site and Loci, in Elements of Evolutionary
 * Genetics, Roberts &amp; Company Publishers, pp. 336-443.
 * </a>
 * </li>
 *
 * </ul>
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public final class LinkageDisequilibrium {

    // Chi-square distribution with 1 degree of freedom - used to compute significance of r2
    private static final ChiSquaredDistribution CHI_SQUARE_1_DF =
            new ChiSquaredDistribution(null, 1);

    // cannot be instantiated
    private LinkageDisequilibrium() {}

    /**
     * Computes the linkage disequilibrium determinant (D) as defined by
     * <a href="http://www.genetics.org/content/49/1/49">Lewontin (1964)</a>.
     *
     * <p>This implementation computes the directional linkage disequilibria (D<SUB>&omega;</SUB>)
     * following formula 1 in
     * <a href="http://www.genetics.org/content/78/3/937.long">Langley &amp; Crow (1974)</a>, and
     * thus requires already polarized frequencies.
     *
     * @param pA  major allele frequency of locus A.
     * @param pB  major allele frequency of locus B.
     * @param pAB haplotype frequency of major alleles carriers (for locus A and B).
     *
     * @return linkage disequilibrium determinant (directional).
     */
    public static double d(final double pA, final double pB, final double pAB) {
        // validate that the frequencies are major
        FrequencyUtils.validateMayorFrequency(pA, () -> "locus A (pA)");
        FrequencyUtils.validateMayorFrequency(pB, () -> "locus A (pB)");
        // validates the frequency range of the haplotype AB
        FrequencyUtils.validateFrequencyRange(pAB);
        return pAB - pA * pB;
    }

    /**
     * Computes the normalized linkage disequilibrium determinant (D') as defined by
     * <a href="http://www.genetics.org/content/49/1/49">Lewontin (1964)</a>.
     *
     * <p>This implementation follows formula 13 in
     * <a href="http://www.sciencedirect.com/science/article/pii/S0040580908000609">
     * VanLiere &amp; Rosenberg (2008)</a>.
     *
     * <p>Note: it requires polarized frequencies (major allele frequency) for computing D using
     * {@link #d(double, double, double)}.
     *
     * @param pA  major allele frequency of locus A.
     * @param pB  major allele frequency of locus B.
     * @param pAB haplotype frequency of major alleles carriers (for locus A and B).
     *
     * @return linkage disequilibrium determinant.
     *
     * @see #d(double, double, double)
     */
    public static double dPrime(final double pA, final double pB, final double pAB) {
        // compute D (also validates frequencies)
        final double d = d(pA, pB, pAB);
        // early termination
        if (d == 0) {
            return d;
        }
        // compute maximum D as described after formula 13 in VanLiere & Rosenberg (2008)
        final double dMax = (d < 0)
                ? FastMath.min(pA * pB, (1 - pA) * (1 - pB))
                : FastMath.min(pA * (1 - pB), (1 - pA) * pB);
        // return D / D' (formula 13)
        return d / dMax;
    }

    /**
     * Computes allele frequency signed correlation (<i>r</i><SUB>&omega;</SUB>).
     *
     * <p>This implementation uses D<SUB>&omega;</SUB> ({@link #d(double, double, double)}) to
     * compute the directionality of the correlation. Thus, it requires polarized frequencies
     * (major allele frequency).
     *
     * @param pA  major allele frequency of locus A.
     * @param pB  major allele frequency of locus B.
     * @param pAB haplotype frequency of major alleles carriers (for locus A and B).
     *
     * @return signed correlation.
     *
     * @see #d(double, double, double)
     */
    public static double rw(final double pA, final double pB, final double pAB) {
        // compute D (and validates frequencies)
        final double d = d(pA, pB, pAB);
        return d / FastMath.sqrt(pA * (1 - pA) * pB * (1 - pB));
    }

    /**
     * Computes allele frequency Pearson's correlation (<i>r</i><SUP>2</SUP>) as defined by
     * <a href="https://link.springer.com/article/10.1007%2FBF01245622">
     * Hill &amp; Robertson (1968)</a>.
     *
     * <p>This implementation follows formula 1 in
     * <a href="http://www.sciencedirect.com/science/article/pii/S0040580908000609">
     * VanLiere &amp; Rosenberg (2008)</a>.
     *
     * <p>Note: it requires polarized frequencies (major allele frequency) for computing D using
     * {@link #d(double, double, double)}.
     *
     * @param pA  major allele frequency of locus A.
     * @param pB  major allele frequency of locus B.
     * @param pAB haplotype frequency of major alleles carriers (for locus A and B).
     *
     * @return pearson correlation.
     */
    public static double r2(final double pA, final double pB, final double pAB) {
        // computes r (and validates frequencies)
        final double absRw = FastMath.abs(rw(pA, pB, pAB));
        // square r
        return absRw * absRw;
    }

    /**
     * Computes allele frequency maximum correlation (<i>r</i><SUP>2</SUP><SUB>max</SUB>) following
     * <a href="http://www.sciencedirect.com/science/article/pii/S0040580908000609">
     * VanLiere &amp; Rosenberg (2008)</a> (formulas 2 and 3).
     *
     * <p>This implementation does not use formula 1 or 4 because they correspond to the space
     * where the frequencies of locus A and B are minor. Thus, it requires polarized frequencies
     * (major allele frequency).
     *
     * @param pA major allele frequency of locus A.
     * @param pB major allele frequency of locus B.
     *
     * @return the maximum correlation.
     */
    public static double r2max(final double pA, final double pB) {
        // validate that the frequencies are major
        FrequencyUtils.validateMayorFrequency(pA, () -> "locus A (pA)");
        FrequencyUtils.validateMayorFrequency(pB, () -> "locus A (pB)");
        if (pA == pB) {
            return 1;
        }
        // polarize locus A and B to simplify the computation using only formula 1
        // see the comment in r2maxPolarized for more information about this logic
        return (pA < pB) ? r2maxPolarized(pA, pB) : r2maxPolarized(pB, pA);
    }

    // helper method for r2max for simplification, computes the maximum correlation with already
    // polarized frequencies based on Table 1 from VanLiere & Rosenberg (2008):
    // 1. We assume that the frequencies are already major, so that reduces the computation to the
    // rows where the columns 'pA < 1/2' and 'pB < 1/2' are 'No' (S2 and S3).
    // 2. Both S2 and S3 are representing the same function, but in one case pA < pB and in the
    // other pB < pA (column 'pA < pB'). Because what is called locus A and locus B is arbitrary,
    // we call locus A the one which satisfy the first condition (pA < pB) and use the formula in
    // S2.
    private static double r2maxPolarized(final double smaller, final double higher) {
        return (smaller * (1 - higher)) / ((1 - smaller) * higher);
    }

    /**
     * Computes allele frequency normalized correlation <i>r</i><SUP>2</SUP>') following
     * <a href="http://www.sciencedirect.com/science/article/pii/S0040580908000609">
     * VanLiere &amp; Rosenberg (2008)</a>
     * (<code><i>r</i><SUP>2</SUP> / <i>r</i><SUP>2</SUP><SUB>max</SUB></code>).
     *
     * <p>Note: it requires polatized frequencies for the computation of r<SUP>2</SUP>
     * ({@link #r2(double, double, double)}) and <i>r</i><SUP>2</SUP><SUB>max</SUB>
     * ({@link #r2max(double, double)}).
     *
     * @param pA  major allele frequency of locus A.
     * @param pB  major allele frequency of locus B.
     * @param pAB haplotype frequency of major alleles carriers (for locus A and B).
     *
     * @return normalized correlation.
     *
     * @see #r2(double, double, double)
     * @see #r2max(double, double)
     */
    public static double r2Prime(final double pA, final double pB, final double pAB) {
        // frequencies are validated in both r2max and r2, so no need to do it here
        final double maxR2 = r2max(pA, pB);
        final double r2 = r2(pA, pB, pAB);
        return r2 / maxR2;
    }


    /**
     * Tests if a value for the allele frequency correlation (<i>r</i><SUP>2</SUP>) is significant
     * based on the number of samples used for its computation
     *
     * <p>This test is based on the property of 2x2 tables described in
     * <a href="https://books.google.com/books/about/Elements_of_Evolutionary_Genetics.html?id=dgNFAQAAIAAJ">
     * Charlesworth &amp; Charlesworth (2012)</a>, formula B8.3.1, and it is based on a Chi-square
     * distribution with 1 degree of freedom.
     *
     * <p>Note that this test cannot be used for small samples (&le; 5).
     *
     * <p>This can be used for test a correlation computed with {@link #r2(double, double, double)}
     * or check if the maximum correlation that can be reached by two locus
     * ({@link #r2max(double, double)}) can be significant.
     *
     * @param r2              correlation to test.
     * @param numberOfSamples number of samples used to compute the correlation (more than 5).
     * @param chiSqrQuantile  the Chi-Squared quantile to asses significance.
     *
     * @return {@code true} if the correlation is significant; {@code false} otherwise.
     *
     * @see #r2SignificantThreshold(int, double)
     */
    public static boolean r2SignificantTest(final double r2, final int numberOfSamples,
            final double chiSqrQuantile) {
        // validate the range of r2
        Verify.validate(r2 >= 0 && r2 <= 1,
                () -> "r2 range should be between 0 and 1: " + numberOfSamples);
        return r2 >= r2SignificantThreshold(numberOfSamples, chiSqrQuantile);
    }

    /**
     * Computes the minimum allele frequency correlation (<i>r</i><SUP>2</SUP>) that could be
     * significant for the provided number of samples used in its computation.
     *
     * <p>This test is based on the property of 2x2 tables described in
     * <a href="https://books.google.com/books/about/Elements_of_Evolutionary_Genetics.html?id=dgNFAQAAIAAJ">
     * Charlesworth &amp; Charlesworth (2012)</a>, formula B8.3.1, and it is based on a Chi-square
     * distribution with 1 degree of freedom.
     *
     * <p>Note that this test cannot be used for small samples (&le; 5).
     *
     * @param numberOfSamples number of samples used to compute the correlation (more than 5).
     * @param chiSqrQuantile  Chi-Squared quantile to asses significance.
     *
     * @return the minimum correlation that can be considered significantly different from
     * {@code 0} for the provided number of samples.
     */
    public static double r2SignificantThreshold(final int numberOfSamples,
            final double chiSqrQuantile) {
        // validate the number of samples
        Verify.validate(numberOfSamples > 5,
                () -> "numberOfSamples should be larger than 1: " + numberOfSamples);
        // quantile is validated in the chi-square implementation
        return CHI_SQUARE_1_DF.inverseCumulativeProbability(chiSqrQuantile) / numberOfSamples;
    }

}
