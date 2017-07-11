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

package org.magicdgs.popgenlib;

import org.openjdk.jmh.annotations.Mode;
import org.openjdk.jmh.results.RunResult;
import org.openjdk.jmh.runner.Runner;
import org.openjdk.jmh.runner.options.Options;
import org.openjdk.jmh.runner.options.OptionsBuilder;
import org.openjdk.jmh.runner.options.TimeValue;
import org.testng.annotations.Test;

import java.util.Collection;
import java.util.Map;
import java.util.concurrent.TimeUnit;
import java.util.stream.Collectors;

/**
 * Base class for all performance tests.
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
@Test(groups = "performance")
public class PopGenLibPerformanceTest extends PopGenLibTest {

    /**
     * Run all the implemented benchmarks for this class, in one thread and one fork, performing GC
     * between each iteration.
     *
     * <p>Note: the benchmark does not use the JVM arguments provided, to override the gradle/IDE
     * defaults.
     *
     * @param mode                  the mode to run the benchmark.
     * @param warmpupIterations     number of warmup iterations (one second per iteration).
     * @param measurementIterations number of iterations for the benchmark code.
     *
     * @return a map with the benchmark name and the score for it to test performance improvements.
     * The time unit returned is microseconds.
     *
     * @throws Exception if the benchmark fails.
     */
    public final Map<String, Double> runBenchmark(final Mode mode,
            final int warmpupIterations, final int measurementIterations) throws Exception {
        final Options opt = new OptionsBuilder()
                .include(this.getClass().getName() + ".*")
                .mode(mode)
                .timeUnit(TimeUnit.MICROSECONDS)
                .warmupTime(TimeValue.seconds(1))
                .warmupIterations(warmpupIterations)
                .measurementTime(TimeValue.seconds(1))
                .measurementIterations(measurementIterations)
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
}
