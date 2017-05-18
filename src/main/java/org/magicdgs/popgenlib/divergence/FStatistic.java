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

import org.magicdgs.popgenlib.diversity.NucleotideDiversity;
import org.magicdgs.popgenlib.utils.Verify;

import java.util.ArrayList;
import java.util.List;

/**
 * Methods to compute divergence using F<sub>ST</sub>
 *
 * In addition, it includes the hierarchical F-statistics: F<sub>RT</sub> and F<sub>RS</sub>.
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public final class FStatistic {

    // cannot be instantiated
    private FStatistic() {}

    /**
     * Computes F<sub>ST</sub> for a pair of populations.
     *
     * @param numberOfSamples1   number of samples (population 1) used for get the frequencies..
     * @param alleleFrequencies1 frequencies for each allele in population 1.
     * @param numberOfSamples1   number of samples (population 2) used for get the frequencies.
     * @param alleleFrequencies2 frequencies for each allele in population 2.
     */
    // TODO: this formula is similar to the PoPoolation2 implementation
    // TODO: but we require the reference and formula in the javadoc
    // TODO: and it should use the general method
    public static double pairwiseFst(final int numberOfSamples1,
            final List<Double> alleleFrequencies1,
            final int numberOfSamples2, final List<Double> alleleFrequencies2) {
        // validate the frequencies
        Verify.nonNull(alleleFrequencies1, () -> "alleleFrequencies1");
        Verify.nonNull(alleleFrequencies2, () -> "alleleFrequencies2");
        Verify.validate(alleleFrequencies1.size() == alleleFrequencies2.size(),
                () -> "different number of alleles the pair of populations");

        // combine the frequencies for the pair of populations
        // TODO: generalize and extract to a method for several populations
        final List<Double> combinedFrequency = new ArrayList<>(alleleFrequencies1.size());
        for (int i = 0; i < alleleFrequencies1.size(); i++) {
            final double averageFrequency =
                    (alleleFrequencies1.get(i) + alleleFrequencies2.get(i)) / 2;
            combinedFrequency.add(averageFrequency);
        }

        // get the nucleotide diversity for the "combined" population
        // combined frequencies represents the total population, and thus its diversity is the
        // pair-wise differences between population
        final double combinedPi = NucleotideDiversity
                .tajimasPi(numberOfSamples1 + numberOfSamples2, combinedFrequency);

        // if it is 0, do not wait time computing the nucleotide diversity by population
        if (combinedPi == 0) {
            return 0;
        }

        // otherwise, we need the average diversity of the two populations
        // TODO: generalize for compute in several populations
        final double pi1 = NucleotideDiversity.tajimasPi(numberOfSamples1, alleleFrequencies1);
        final double pi2 = NucleotideDiversity.tajimasPi(numberOfSamples2, alleleFrequencies1);
        final double average = (pi1 + pi2) / 2;

        // now the diversity between-within population percentage defines Fst
        return (combinedPi - average) / combinedPi;

    }


    // TODO: documentation
    public static double fst(final List<Integer> numberOfSamplesPerPopulation,
            final List<List<Double>> frequenciesPerPopulation) {
        // validate lists
        Verify.nonEmpty(numberOfSamplesPerPopulation, () -> "numberOfSamplesPerPopulation");
        Verify.nonEmpty(frequenciesPerPopulation, () -> "frequenciesPerPopulation");

        final int nPopulations = numberOfSamplesPerPopulation.size();
        // validate that each sample have a value
        Verify.validate(nPopulations == frequenciesPerPopulation.size(),
                () -> "the number of values in both lists should be equal to the number of populations");

        final int nAlleles = frequenciesPerPopulation.get(0).size();
        // validates that all the populations have the
        frequenciesPerPopulation.forEach(freqs -> Verify.validate(nAlleles == freqs.size(),
                () -> "different number of alleles per population were found"));

        // combine the frequencies for all the populations
        final List<Double> combindedFrequencies = new ArrayList<>(nAlleles);
        for (int i = 0; i < nAlleles; i++) {
            // it should be final for use inside the stream
            final int index = i;
            final double freqSum = frequenciesPerPopulation.stream()
                    .mapToDouble(freqs -> freqs.get(index)).sum();
            combindedFrequencies.add(freqSum / nPopulations);
        }

        // get the nucleotide diversity for the "combined" population
        // combined frequencies represents the total population, and thus its diversity is the
        // pair-wise differences between population
        final int totalSumples = numberOfSamplesPerPopulation.stream().reduce(0, Integer::sum);
        final double betweenPi = NucleotideDiversity.tajimasPi(totalSumples, combindedFrequencies);

        // if it is 0, do not wait time computing the nucleotide diversity by population
        if (betweenPi == 0) {
            return 0;
        }

        // otherwise, we need the average diversity of the two populations
        double averageWithinPi = 0;
        for (int i = 0; i < nAlleles; i++) {
            averageWithinPi += NucleotideDiversity.tajimasPi(
                    numberOfSamplesPerPopulation.get(i),
                    frequenciesPerPopulation.get(i));
        }
        averageWithinPi /= nPopulations;

        // formula 3 in Hudson et al. 1992.
        return 1 - (averageWithinPi / betweenPi);
    }

}
