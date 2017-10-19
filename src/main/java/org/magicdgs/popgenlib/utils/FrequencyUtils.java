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

import org.apache.commons.math3.util.Pair;

import java.math.BigDecimal;
import java.math.MathContext;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.function.Supplier;
import java.util.stream.Collectors;

/**
 * Utility methods for working with frequencies. If frequencies are not valid, these methods throw
 * an {@link IllegalFrequencyException}.
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class FrequencyUtils {

    /** Maximum frequency. */
    public static final Double FREQUENCY_ONE = 1d;

    /** Minimum frequency. */
    public static final Double FREQUENCY_ZERO = 0d;

    // minor frequency trheshold for checking non-minor alleles
    // minor allele should be always lower or equals than 0.5
    private static final double MINOR_FREQUENCY_THRESHOLD = 0.5d;

    // cannot be instantiated
    private FrequencyUtils() {}

    /**
     * Validates that the frequencies are correct:
     * <ul>
     * <li>Frequency list should not be empty.</li>
     * <li>Frequencies should be in the range [{@link #FREQUENCY_ZERO},
     * {@link #FREQUENCY_ONE}].</li>
     * <li>Frequencies should sum exactly {@link #FREQUENCY_ONE}.</li>
     * </ul>
     *
     * @param frequencies frequencies to validate.
     *
     * @throws IllegalFrequencyException if frequencies are not valid.
     */
    public static void validateFrequencies(final List<Double> frequencies) {
        try {
            final double sum = Verify.nonEmpty(frequencies, () -> "frequencies")
                    .stream().mapToDouble(FrequencyUtils::validateFrequencyRange)
                    .sum();
            Verify.validate(FREQUENCY_ONE.equals(sum),
                    () -> String.format("Frequencies should sum 1 but found %s: %s",
                            sum, frequencies));
        } catch (IllegalArgumentException e) {
            throw new IllegalFrequencyException(e);
        }
    }

    /**
     * Validates the range of the frequency ([{@link #FREQUENCY_ZERO}, {@link #FREQUENCY_ONE}]).
     *
     * <p>Note: use {@link #validateFrequencies(List)} for checking if all the frequencies are
     * within range and sum exactly {@link #FREQUENCY_ONE}.
     *
     * @param freq the frequency to validate.
     *
     * @return the same frequency.
     *
     * @throws IllegalFrequencyException if the frequency is not within the range.
     */
    public static double validateFrequencyRange(final Double freq) {
        try {
            Verify.validate(freq >= FREQUENCY_ZERO && freq <= FREQUENCY_ONE,
                    () -> String.format("Frequencies out of range [%s, %s]: %s",
                            FREQUENCY_ZERO, FREQUENCY_ONE, freq));
            return freq;
        } catch (IllegalArgumentException e) {
            throw new IllegalFrequencyException(e);
        }
    }

    /**
     * Validates that the range of the frequency ({@link #validateFrequencyRange(Double)} and
     * that the frequencies corresponds to a minor frequency for a binary trait
     * ({@code freq >= 1/2}).
     *
     * @param freq the frequency to validate.
     * @param msg message to add to the error message for non-minor frequencies.
     *
     * @return the same frequency.
     *
     * @throws IllegalFrequencyException if the frequency is not within the range or minor.
     */
    public static double validateMayorFrequency(final Double freq, final Supplier<String> msg) {
        validateFrequencyRange(freq);
        try {
            Verify.validate(freq >= MINOR_FREQUENCY_THRESHOLD, () -> "Non-minor: " + msg.get());
            return freq;
        } catch (IllegalArgumentException e) {
            throw new IllegalFrequencyException(e);
        }
    }

    /**
     * Gets total counts and frequencies from a list of counts.
     *
     * <p>Note: Original counts could be recovered by multiplying each of the frequencies by the
     * total number of counts.
     *
     * @param counts list of counts for each category.
     *
     * @return total counts and frequencies.
     */
    public static Pair<Integer, List<Double>> countsToFrequencies(final List<Integer> counts) {
        try {
            // verification and counting the total
            Verify.nonEmpty(counts, () -> "allele counts");
            final int total = counts.stream()
                    .mapToInt(i -> {
                        Verify.validate(i >= 0, () -> "counts should be larger than 0: " + i);
                        return i;
                    }).sum();
            Verify.validate(total > 0, () -> "all counts are zero");

            // compute the frequencies using BigDecimal operations
            final BigDecimal totalBig = new BigDecimal(total);
            final List<Double> freqs = new ArrayList<>(counts.size());
            for (final double c : counts) {
                if (c == 0) {
                    freqs.add(FREQUENCY_ZERO);
                } else {
                    freqs.add(new BigDecimal(c).divide(totalBig, MathContext.DECIMAL64)
                            .doubleValue());
                }
            }

            // return a pair
            return Pair.create(total, freqs);
        } catch (IllegalArgumentException e) {
            throw new IllegalFrequencyException(e);
        }
    }

    /**
     * Sort the frequencies from major to minor.
     *
     * @param frequencies the frequencies to sort.
     *
     * @return a new list with the sorted frequencies.
     *
     * @throws IllegalFrequencyException if frequencies are invalid.
     */
    public static List<Double> sortFrequencies(final List<Double> frequencies) {
        validateFrequencies(frequencies);
        return frequencies.stream().sorted(Collections.reverseOrder())
                .collect(Collectors.toList());
    }

    /** Exception throw for invalid frequencies. */
    public static final class IllegalFrequencyException extends IllegalArgumentException {

        /** Wraps an {@link IllegalArgumentException} into a frequency exception. */
        public IllegalFrequencyException(final IllegalArgumentException exception) {
            super(exception.getMessage());
        }
    }
}
