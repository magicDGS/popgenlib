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

import java.util.Collection;
import java.util.function.Supplier;

/**
 * Various utils for validation of parameters.
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public final class Verify {

    // cannot be instantiated
    private Verify() {}

    /**
     * Verifies if the condition is {@code true}.
     *
     * @param condition condition to evaluate.
     * @param msg       message supplier in case of exception.
     *
     * @throws IllegalArgumentException if the condition is {@code false}.
     */
    public static void validate(final boolean condition, final Supplier<String> msg) {
        if (!condition) {
            throw new IllegalArgumentException(msg.get());
        }
    }

    /**
     * Verifies if the object is non-null.
     *
     * @param object object to verify.
     * @param msg    message supplier in case of exception.
     *
     * @return the same object if it pass the verification.
     *
     * @throws IllegalArgumentException if the object is {@code null}.
     */
    public static <T> T nonNull(final T object, final Supplier<String> msg) {
        validate(object != null, msg);
        return object;
    }

    /**
     * Verifies if the collection is not {@code null}, not empty and their contents are not
     * {@code null}.
     *
     * @param collection collection to verify.
     * @param msg        message supplier in case of exception.
     *
     * @return the same object if it pass the verification.
     *
     * @throws IllegalArgumentException if the collection is {@code null}, empty or have
     *                                  {@code null} values.
     */
    public static <T> Collection<T> nonEmpty(final Collection<T> collection,
            final Supplier<String> msg) {
        nonNull(collection, () -> "null collection: " + msg.get());
        validate(!collection.isEmpty(), () -> "empty collection: " + msg.get());
        collection.forEach(value -> nonNull(value, () -> "null value: " + msg.get()));
        return collection;
    }
}
