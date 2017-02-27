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

import org.magicdgs.popgenlib.utils.Verify;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.cache.CacheBuilder;
import com.google.common.cache.CacheLoader;
import com.google.common.cache.LoadingCache;
import com.google.common.util.concurrent.UncheckedExecutionException;

import java.util.stream.IntStream;

/**
 * Class for efficient computation of Watterson's &theta;. This class maintain cached values for
 * the formula denominator, being more efficient in the case that the same number of sample is
 * expected for each calculation.
 *
 * <p>This class is more efficient that using the non-cached variant for cases
 * where calls are expected to do not vary in the number of samples too much. As an example, if
 * &theta;<SUB>W</SUB> should be computed in several windows for the same number of samples,
 * this will avoid recomputing the denominator for each of them.
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 * @see NucleotideDiversity#wattersonsTheta(int, int)
 */
public final class WattersonsThetaCalculator {

    // denominator loader for get the loading cache
    final static CacheLoader<Integer, Double> DENOMINATOR_LOADER =
            new CacheLoader<Integer, Double>() {
                @Override
                public Double load(final Integer numberOfSamples) {
                    Verify.validate(numberOfSamples > 1,
                            () -> "numberOfSamples should be at least 2: " + numberOfSamples);
                    return NucleotideDiversity.wattersonsDenominatorApproximation(numberOfSamples);
                }
            };

    // this is the cache of values for the denominator
    @VisibleForTesting
    final LoadingCache<Integer, Double> denominatorCache;

    /**
     * Creates an instance with default caching.
     *
     * @see CacheBuilder#newBuilder() for default values caching parameters.
     */
    public WattersonsThetaCalculator() {
        this(CacheBuilder.newBuilder());
    }

    /**
     * Creates an instance using the cache provided by the builder.
     *
     * @param cacheBuilder builder for the cached.
     */
    public WattersonsThetaCalculator(final CacheBuilder<Object, Object> cacheBuilder) {
        this.denominatorCache = cacheBuilder.build(DENOMINATOR_LOADER);
    }

    /**
     * Computes Watterson's &theta; using cached denominator.
     *
     * @param numberOfSamples          number of samples used to get the number of segregating
     *                                 sites.
     * @param numberOfSegregatingSites number of segregating sites in the sample.
     * @see NucleotideDiversity#wattersonsTheta(int, int)
     */
    public double wattersonsTheta(final int numberOfSamples,
            final int numberOfSegregatingSites) {
        Verify.validate(numberOfSegregatingSites >= 0,
                () -> "Number of segregating sites should be 0 or a positive integer: "
                        + numberOfSegregatingSites);
        try {
            // use the cached a1 for Watterson's theta
            return numberOfSegregatingSites / denominatorCache.getUnchecked(numberOfSamples);
        } catch (final UncheckedExecutionException e) {
            // the cache throw an unchecked exception, this should be an IllegalArgumentException
            // because no other unchecked exceptions are expected
            // if so, we should handle them properly in NucleotideDiversity
            throw (IllegalArgumentException) e.getCause();
        }
    }
}
