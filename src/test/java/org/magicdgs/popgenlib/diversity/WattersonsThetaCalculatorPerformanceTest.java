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

import org.magicdgs.popgenlib.PopGenLibPerformanceTest;

import org.openjdk.jmh.annotations.Benchmark;
import org.openjdk.jmh.annotations.Mode;
import org.openjdk.jmh.annotations.Scope;
import org.openjdk.jmh.annotations.State;
import org.openjdk.jmh.infra.Blackhole;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Map;

/**
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class WattersonsThetaCalculatorPerformanceTest extends PopGenLibPerformanceTest {

    @Test(singleThreaded = true)
    public void testFasterCachedVersion() throws Exception {
        final Map<String, Double> avgTime = runBenchmark(Mode.AverageTime, 2, 20);
        final double cached = avgTime.get(this.getClass().getName() + ".MicroBenchmark.cached");
        final double nonCached =
                avgTime.get(this.getClass().getName() + ".MicroBenchmark.nonCached");
        // we expect improvement in the performance, even if it is small
        if (cached > nonCached) {
            Assert.fail(String.format(
                    "No improvement in performance: %s ms/op (cached) vs. %s ms/op (non cached)",
                    cached, nonCached
            ));
        }
    }

    // TODO: maybe extract this benchmarks to a different module and run performance tests there
    @State(Scope.Thread)
    public static class MicroBenchmark {

        @Benchmark
        public void nonCached(final MicroBenchmark state, final Blackhole bh) {
            for (int i = 0; i < 100000; i++) {
                bh.consume(NucleotideDiversity.wattersonsTheta(1000, 2));
            }
        }

        @Benchmark
        public void cached(final MicroBenchmark state, final Blackhole bh) {
            final WattersonsThetaCalculator calculator = new WattersonsThetaCalculator();
            for (int i = 0; i < 100000; i++) {
                bh.consume(calculator.wattersonsTheta(1000, 2));
            }
        }
    }
}
