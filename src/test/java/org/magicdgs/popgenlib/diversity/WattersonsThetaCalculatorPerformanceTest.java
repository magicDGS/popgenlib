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

import org.magicdgs.popgenlib.PopGenLibTest;

import org.openjdk.jmh.annotations.Benchmark;
import org.openjdk.jmh.annotations.Mode;
import org.openjdk.jmh.annotations.Scope;
import org.openjdk.jmh.annotations.State;
import org.openjdk.jmh.infra.Blackhole;
import org.openjdk.jmh.results.RunResult;
import org.openjdk.jmh.runner.Runner;
import org.openjdk.jmh.runner.options.Options;
import org.openjdk.jmh.runner.options.OptionsBuilder;
import org.openjdk.jmh.runner.options.TimeValue;
import org.openjdk.jmh.runner.options.VerboseMode;
import org.testng.Assert;
import org.testng.annotations.AfterSuite;
import org.testng.annotations.Test;

import java.util.Collection;
import java.util.Map;
import java.util.concurrent.TimeUnit;
import java.util.stream.Collectors;

/**
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class WattersonsThetaCalculatorPerformanceTest extends PopGenLibTest {

    // TODO: move to PopGenLibTest to run all benchmarks??
    public final Map<String, Double> runBenchmark(final Mode mode) throws Exception {
        final Options opt = new OptionsBuilder()
                .include(this.getClass().getName() + ".*")
                .mode(mode)
                .timeUnit(TimeUnit.MICROSECONDS)
                .warmupTime(TimeValue.seconds(1))
                .warmupIterations(2)
                .measurementTime(TimeValue.seconds(1))
                .measurementIterations(10)
                .threads(1)
                .forks(1)
                .shouldFailOnError(true)
                .shouldDoGC(true)
                // do not run with JVM to override the gradle/IDE defaults
                .jvmArgs()
                .build();

        // run the benchmark
        final Collection<RunResult> results = new Runner(opt).run();
        // creates a map with the benchmark name and the score
        return results.stream().collect(Collectors.toMap(
                result -> result.getParams().getBenchmark(),
                result -> result.getPrimaryResult().getScore()));
    }

    @Test(singleThreaded = true)
    public void testSpeedComparedToNonCached() throws Exception {
        final Map<String, Double> avgTime = runBenchmark(Mode.AverageTime);
        final double cached = avgTime.get(this.getClass().getName() + ".MicroBenchmark.cached");
        final double nonCached = avgTime.get(this.getClass().getName() + ".MicroBenchmark.nonCached");
        // we expect improvement in the performance, even if it is small
        if (cached > nonCached) {
            Assert.fail(String.format("No improvement in performance: %s ms/op (cached) vs. %s ms/op (non cached)",
                cached, nonCached
            ));
        }
    }

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
