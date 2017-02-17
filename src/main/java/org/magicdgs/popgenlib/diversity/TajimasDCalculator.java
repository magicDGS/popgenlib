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

import com.google.common.cache.CacheBuilder;
import com.google.common.cache.CacheLoader;
import com.google.common.cache.LoadingCache;
import org.apache.commons.math3.util.Pair;

/**
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
// TODO: document
public class TajimasDCalculator {

    // this is the cache for a1
    private final LoadingCache<Integer, Double> a1cache;
    // this is the cache for the constants e1 and e2 (includes a1 cache)
    private final LoadingCache<Integer, Pair<Double, Double>> constantCache;


    /**
     * Creates an instance using the cache provided by the builder.
     *
     * @param cacheBuilder builder for the cached.
     */
    public TajimasDCalculator(final CacheBuilder<Object, Object> cacheBuilder) {
        this.a1cache = cacheBuilder.build(WattersonsThetaCalculator.DENOMINATOR_LOADER);
        // constant cache use the a1 cache to avoid recomputation
        this.constantCache = cacheBuilder.build(new CacheLoader<Integer, Pair<Double, Double>>() {
            @Override
            public Pair<Double, Double> load(Integer key) throws Exception {
                return TajimasD.varianceConstants(a1cache.getUnchecked(key), key);
            }
        });
    }


    public double tajimasD(final int numberOfSamples, final int piEstimate, final int segregatingSites) {
        return TajimasD.tajimasD(numberOfSamples, piEstimate,
                () -> wattersonsTheta(numberOfSamples, segregatingSites),
                constantCache.getUnchecked(numberOfSamples));
    }

    /**
     * Computes Watterson's &theta; using cached denominator.
     *
     * @param numberOfSamples          number of samples used to get the number of segregating
     *                                 sites.
     * @param numberOfSegregatingSites number of segregating sites in the sample.
     *
     * @see NucleotideDiversity#wattersonsTheta(int, int)
     */
    public double wattersonsTheta(final int numberOfSamples,
            final int numberOfSegregatingSites) {
        Verify.validate(numberOfSamples > 1,
                () -> "numberOfSamples should be at least 2 for computing Watterson's Theta: "
                        + numberOfSamples);

        // - S is the number of segregating sites
        // - a is the denominator of Watterson's theta (using the cache)
        // we use the unchecked cache because the denominator loader does not throws any checked exception
        return NucleotideDiversity.wattersonsTheta(numberOfSegregatingSites,
                a1cache.getUnchecked(numberOfSamples));
    }

}
