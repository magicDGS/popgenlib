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

import org.magicdgs.popgenlib.PopGenLibTest;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.function.Supplier;

/**
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class VerifyUnitTest extends PopGenLibTest {

    private final static Supplier<String> EMPTY_SUPPLIER = () -> "";

    @Test
    public void testValidate() throws Exception {
        Assert.assertThrows(IllegalArgumentException.class,
                () -> Verify.validate(false, EMPTY_SUPPLIER));
        Verify.validate(true, EMPTY_SUPPLIER);
    }

    @Test
    public void testNonNull() throws Exception {
        // null throws
        Assert.assertThrows(IllegalArgumentException.class,
                () -> Verify.nonNull(null, EMPTY_SUPPLIER));
        // integer 1 do not throws
        final Integer integer = 1;
        Assert.assertSame(integer, Verify.nonNull(integer, EMPTY_SUPPLIER));
    }

    @DataProvider(name = "emptyCollections")
    public Object[][] getInvalidCollections() {
        return new Object[][] {
                {null},
                {Collections.emptySet()},
                {Collections.singleton(null)},
                {Arrays.asList(1, null)},
                {Arrays.asList(null, 1)}
        };
    }


    @Test(dataProvider = "emptyCollections", expectedExceptions = IllegalArgumentException.class)
    public void testInvalidCollections(final Collection<?> collection) throws Exception {
        Verify.nonEmpty(collection, EMPTY_SUPPLIER);
    }

    @DataProvider(name = "nonEmptyCollections")
    public Object[][] getValidCollections() {
        return new Object[][] {
                {Collections.singleton(1)},
                {Collections.singletonList(1)},
                {Arrays.asList(1, 2, 3)}
        };
    }

    @Test(dataProvider = "nonEmptyCollections")
    public void testNonEmpty(final Collection<?> collection) throws Exception {
        Assert.assertSame(Verify.nonEmpty(collection, EMPTY_SUPPLIER), collection);
    }
}