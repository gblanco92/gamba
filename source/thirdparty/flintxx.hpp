/* Copyright 2001-2003, 2006, 2008, 2011-2015, 2018 Free Software
Foundation, Inc.

This file was modified from the GNU MP Library.

The GNU MP Library is free software; you can redistribute it and/or modify
it under the terms of either:

  * the GNU Lesser General Public License as published by the Free
    Software Foundation; either version 3 of the License, or (at your
    option) any later version.

or

  * the GNU General Public License as published by the Free Software
    Foundation; either version 2 of the License, or (at your option) any
    later version.

or both in parallel, as here.

The GNU MP Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received copies of the GNU General Public License and the
GNU Lesser General Public License along with the GNU MP Library.  If not,
see https://www.gnu.org/licenses/.  */

#pragma once

#include <stdexcept>
#include <string>

// clang-format off
#include <flint/flint.h>
#include <gmp.h>
#include <flint/fmpq.h>
#include <flint/fmpz.h>
// clang-format on

#ifdef __clang__
#pragma clang diagnostic ignored "-Wdouble-promotion"
#elif __GNUC__
#pragma GCC diagnostic ignored "-Wdouble-promotion"
#endif

// NOLINTBEGIN(google-runtime-int)
// NOLINTBEGIN(bugprone-macro-parentheses)
// NOLINTBEGIN(google-explicit-constructor,hicpp-explicit-conversions)
// NOLINTBEGIN(cert-dcl51-cpp,cert-dcl37-c,bugprone-reserved-identifier)

#define __FLINTXX_CONSTANT(X)      __builtin_constant_p(X)
#define __FLINTXX_CONSTANT_TRUE(X) (__FLINTXX_CONSTANT(X) && (X))

/**************** Function objects ****************/
/* Any evaluation of a __flint_expr ends up calling one of these functions
   all intermediate functions being inline, the evaluation should optimize
   to a direct call to the relevant function, thus yielding no overhead
   over the C interface. */

struct __flint_unary_plus
{
    static void eval(fmpz* z, fmpz const* w) { fmpz_set(z, w); }

    static void eval(fmpq* q, fmpq const* r) { fmpq_set(q, r); }
};

struct __flint_unary_minus
{
    static void eval(fmpz* z, fmpz const* w) { fmpz_neg(z, w); }

    static void eval(fmpq* q, fmpq const* r) { fmpq_neg(q, r); }
};

struct __flint_unary_com
{
    static void eval(fmpz* z, fmpz const* w) { fmpz_complement(z, w); }
};

struct __flint_binary_plus
{
    static void eval(fmpz* z, fmpz const* w, fmpz const* v)
    {
        fmpz_add(z, w, v);
    }

    static void eval(fmpz* z, fmpz const* w, unsigned long int l)
    {
        // Ideally, those checks should happen earlier so that the tree
        // generated for a+0+b would just be sum(a,b).
        if (__FLINTXX_CONSTANT(l) && l == 0)
            fmpz_set(z, w);
        else
            fmpz_add_ui(z, w, l);
    }

    static void eval(fmpz* z, unsigned long int l, fmpz const* w)
    {
        eval(z, w, l);
    }

    static void eval(fmpz* z, fmpz const* w, signed long int l)
    {
        if (l >= 0)
            eval(z, w, static_cast<unsigned long>(l));
        else
            fmpz_sub_ui(z, w, -static_cast<unsigned long>(l));
    }

    static void eval(fmpz* z, signed long int l, fmpz const* w)
    {
        eval(z, w, l);
    }

    static void eval(fmpz* z, fmpz const* w, double d)
    {
        fmpz_add_si(z, w, static_cast<slong>(d));
    }

    static void eval(fmpz* z, double d, fmpz const* w) { eval(z, w, d); }

    static void eval(fmpq* q, fmpq const* r, fmpq const* s)
    {
        fmpq_add(q, r, s);
    }

    static void eval(fmpq* q, fmpq const* r, unsigned long int l)
    {
        if (__FLINTXX_CONSTANT(l) && l == 0)
        {
            fmpq_set(q, r);
        }
        else if (__FLINTXX_CONSTANT(l) && l == 1)
        {
            fmpz_add(fmpq_numref(q), fmpq_numref(r), fmpq_denref(r));

            fmpz_set(fmpq_denref(q), fmpq_denref(r));
        }
        else
        {
            if (q == r)
            {
                fmpz_addmul_ui(fmpq_numref(q), fmpq_denref(q), l);
            }
            else
            {
                fmpz_mul_ui(fmpq_numref(q), fmpq_denref(r), l);
                fmpz_add(fmpq_numref(q), fmpq_numref(q), fmpq_numref(r));
                fmpz_set(fmpq_denref(q), fmpq_denref(r));
            }
        }
    }

    static void eval(fmpq* q, unsigned long int l, fmpq const* r)
    {
        eval(q, r, l);
    }

    static inline void eval(fmpq* q, fmpq const* r, signed long int l);

    // defined after __flint_binary_minus
    static void eval(fmpq* q, signed long int l, fmpq const* r)
    {
        eval(q, r, l);
    }

    static void eval(fmpq* q, fmpq const* r, double d)
    {
        fmpq_add_si(q, r, static_cast<slong>(d));
    }

    static void eval(fmpq* q, double d, fmpq const* r) { eval(q, r, d); }

    static void eval(fmpq* q, fmpq const* r, fmpz const* z)
    {
        if (q == r)
        {
            fmpz_addmul(fmpq_numref(q), fmpq_denref(q), z);
        }
        else
        {
            fmpz_mul(fmpq_numref(q), fmpq_denref(r), z);
            fmpz_add(fmpq_numref(q), fmpq_numref(q), fmpq_numref(r));
            fmpz_set(fmpq_denref(q), fmpq_denref(r));
        }
    }

    static void eval(fmpq* q, fmpz const* z, fmpq const* r) { eval(q, r, z); }
};

struct __flint_binary_minus
{
    static void eval(fmpz* z, fmpz const* w, fmpz const* v)
    {
        fmpz_sub(z, w, v);
    }

    static void eval(fmpz* z, fmpz const* w, unsigned long int l)
    {
        if (__FLINTXX_CONSTANT(l) && l == 0)
        {
            fmpz_set(z, w);
        }
        else
        {
            fmpz_sub_ui(z, w, l);
        }
    }

    static void eval(fmpz* z, unsigned long int l, fmpz const* w)
    {
        if (__FLINTXX_CONSTANT(l) && l == 0)
        {
            fmpz_neg(z, w);
        }
        else
        {
            fmpz_set(z, w);
            fmpz_neg(z, z);
            eval(z, w, l);
        }
    }

    static void eval(fmpz* z, fmpz const* w, signed long int l)
    {
        if (l >= 0)
            eval(z, w, static_cast<unsigned long>(l));
        else
            fmpz_add_ui(z, w, -static_cast<unsigned long>(l));
    }

    static void eval(fmpz* z, signed long int l, fmpz const* w)
    {
        if (l >= 0)
        {
            eval(z, static_cast<unsigned long>(l), w);
        }
        else
        {
            fmpz_add_ui(z, w, -static_cast<unsigned long>(l));
            fmpz_neg(z, z);
        }
    }

    static void eval(fmpz* z, fmpz const* w, double d)
    {
        fmpz_sub_si(z, w, static_cast<slong>(d));
    }

    static void eval(fmpz* z, double d, fmpz const* w)
    {
        eval(z, static_cast<slong>(d), w);
    }

    static void eval(fmpq* q, fmpq const* r, fmpq const* s)
    {
        fmpq_sub(q, r, s);
    }

    static void eval(fmpq* q, fmpq const* r, unsigned long int l)
    {
        if (__FLINTXX_CONSTANT(l) && l == 0)
        {
            fmpq_set(q, r);
        }
        else if (__FLINTXX_CONSTANT(l) && l == 1)
        {
            fmpz_sub(fmpq_numref(q), fmpq_numref(r), fmpq_denref(r));

            fmpz_set(fmpq_denref(q), fmpq_denref(r));
        }
        else
        {
            if (q == r)
            {
                fmpz_submul_ui(fmpq_numref(q), fmpq_denref(q), l);
            }
            else
            {
                fmpz_mul_ui(fmpq_numref(q), fmpq_denref(r), l);
                fmpz_sub(fmpq_numref(q), fmpq_numref(r), fmpq_numref(q));
                fmpz_set(fmpq_denref(q), fmpq_denref(r));
            }
        }
    }

    static void eval(fmpq* q, unsigned long int l, fmpq const* r)
    {
        eval(q, r, l);
        fmpq_neg(q, q);
    }

    static void eval(fmpq* q, fmpq const* r, signed long int l)
    {
        if (l >= 0)
            eval(q, r, static_cast<unsigned long>(l));
        else
            __flint_binary_plus::eval(q, r, -static_cast<unsigned long>(l));
    }

    static void eval(fmpq* q, signed long int l, fmpq const* r)
    {
        eval(q, r, l);
        fmpq_neg(q, q);
    }

    static void eval(fmpq* q, fmpq const* r, double d)
    {
        fmpq_sub_si(q, r, static_cast<slong>(d));
    }

    static void eval(fmpq* q, double d, fmpq const* r)
    {
        eval(q, r, -d);
        fmpq_neg(q, q);
    }

    static void eval(fmpq* q, fmpq const* r, fmpz const* z)
    {
        if (q == r)
        {
            fmpz_submul(fmpq_numref(q), fmpq_denref(q), z);
        }
        else
        {
            fmpz_mul(fmpq_numref(q), fmpq_denref(r), z);
            fmpz_sub(fmpq_numref(q), fmpq_numref(r), fmpq_numref(q));
            fmpz_set(fmpq_denref(q), fmpq_denref(r));
        }
    }

    static void eval(fmpq* q, fmpz const* z, fmpq const* r)
    {
        eval(q, r, z);
        fmpq_neg(q, q);
    }
};

// defined here so it can reference __flint_binary_minus
inline void __flint_binary_plus::eval(fmpq* q, fmpq const* r, signed long int l)
{
    if (l >= 0)
        eval(q, r, static_cast<unsigned long>(l));
    else
        __flint_binary_minus::eval(q, r, -static_cast<unsigned long>(l));
}

struct __flint_binary_lshift
{
    static void eval(fmpz* z, fmpz const* w, mp_bitcnt_t l)
    {
        if (__FLINTXX_CONSTANT(l) && (l == 0))
            fmpz_set(z, w);
        else
            fmpz_mul_2exp(z, w, l);
    }

    static void eval(fmpq* q, fmpq const* r, mp_bitcnt_t l)
    {
        if (__FLINTXX_CONSTANT(l) && (l == 0))
            fmpq_set(q, r);
        else
            fmpq_mul_2exp(q, r, l);
    }
};

struct __flint_binary_rshift
{
    static void eval(fmpz* z, fmpz const* w, mp_bitcnt_t l)
    {
        if (__FLINTXX_CONSTANT(l) && (l == 0))
            fmpz_set(z, w);
        else
            fmpz_fdiv_q_2exp(z, w, l);
    }

    static void eval(fmpq* q, fmpq const* r, mp_bitcnt_t l)
    {
        if (__FLINTXX_CONSTANT(l) && (l == 0))
            fmpq_set(q, r);
        else
            fmpq_div_2exp(q, r, l);
    }
};

struct __flint_binary_multiplies
{
    static void eval(fmpz* z, fmpz const* w, fmpz const* v)
    {
        fmpz_mul(z, w, v);
    }

    static void eval(fmpz* z, fmpz const* w, unsigned long int l)
    {
        fmpz_mul_ui(z, w, l);
    }

    static void eval(fmpz* z, unsigned long int l, fmpz const* w)
    {
        eval(z, w, l);
    }

    static void eval(fmpz* z, fmpz const* w, signed long int l)
    {
        if (__FLINTXX_CONSTANT_TRUE(l >= 0))
        {
            eval(z, w, static_cast<unsigned long>(l));
        }
        else if (__FLINTXX_CONSTANT_TRUE(l <= 0))
        {
            eval(z, w, -static_cast<unsigned long>(l));
            fmpz_neg(z, z);
        }
        else
        {
            fmpz_mul_si(z, w, l);
        }
    }

    static void eval(fmpz* z, signed long int l, fmpz const* w)
    {
        eval(z, w, l);
    }

    static void eval(fmpz* z, fmpz const* w, double d)
    {
        fmpz_mul_si(z, w, static_cast<slong>(d));
    }

    static void eval(fmpz* z, double d, fmpz const* w) { eval(z, w, d); }

    static void eval(fmpq* q, fmpq const* r, fmpq const* s)
    {
        fmpq_mul(q, r, s);
    }

    static void eval(fmpq* q, fmpq const* r, unsigned long int l)
    {
        fmpq_mul_ui(q, r, l);
    }

    static void eval(fmpq* q, unsigned long int l, fmpq const* r)
    {
        eval(q, r, l);
    }

    static void eval(fmpq* q, fmpq const* r, signed long int l)
    {
        if (__FLINTXX_CONSTANT_TRUE(l >= 0))
        {
            eval(q, r, static_cast<unsigned long>(l));
        }
        else if (__FLINTXX_CONSTANT_TRUE(l <= 0))
        {
            eval(q, r, -static_cast<unsigned long>(l));
            fmpq_neg(q, q);
        }
        else
        {
            fmpq_mul_si(q, r, l);
        }
    }

    static void eval(fmpq* q, signed long int l, fmpq const* r)
    {
        eval(q, r, l);
    }

    static void eval(fmpq* q, fmpq const* r, double d)
    {
        fmpq_mul_si(q, r, static_cast<slong>(d));
    }

    static void eval(fmpq* q, double d, fmpq const* r) { eval(q, r, d); }
};

struct __flint_binary_divides
{
    static void eval(fmpz* z, fmpz const* w, fmpz const* v)
    {
        fmpz_tdiv_q(z, w, v);
    }

    static void eval(fmpz* z, fmpz const* w, unsigned long int l)
    {
        fmpz_tdiv_q_ui(z, w, l);
    }

    static void eval(fmpz* z, unsigned long int l, fmpz const* w)
    {
        if (fmpz_sgn(w) >= 0)
        {
            if (fmpz_abs_fits_ui(w))
                fmpz_set_ui(z, l / fmpz_get_ui(w));
            else
                fmpz_set_ui(z, 0);
        }
        else
        {
            fmpz_neg(z, w);

            if (fmpz_abs_fits_ui(z))
            {
                fmpz_set_ui(z, l / fmpz_get_ui(z));
                fmpz_neg(z, z);
            }
            else
            {
                fmpz_set_ui(z, 0);
            }
        }
    }

    static void eval(fmpz* z, fmpz const* w, signed long int l)
    {
        if (l >= 0)
        {
            eval(z, w, static_cast<unsigned long>(l));
        }
        else
        {
            eval(z, w, -static_cast<unsigned long>(l));
            fmpz_neg(z, z);
        }
    }

    static void eval(fmpz* z, signed long int l, fmpz const* w)
    {
        if (fmpz_abs_fits_ui(w))
        {
            fmpz_set_si(z, l / fmpz_get_si(w));
        }
        else
        {
            fmpz_t tmp;
            fmpz_set_si(tmp, l);
            /* if w is bigger than a long then the quotient must be zero,
               unless l==LONG_MIN and w==-LONG_MIN in which case the
               quotient is -1 */
            fmpz_set_si(z, (fmpz_cmpabs(w, tmp) == 0 ? -1 : 0));
            fmpz_clear(tmp);
        }
    }

    static void eval(fmpz* z, fmpz const* w, double d)
    {
        __flint_binary_multiplies::eval(z, w, static_cast<slong>(1.0 / d));
    }

    static void eval(fmpz* z, double d, fmpz const* w)
    {
        __flint_binary_multiplies::eval(z, static_cast<slong>(1.0 / d), w);
    }

    static void eval(fmpq* q, fmpq const* r, fmpq const* s)
    {
        fmpq_div(q, r, s);
    }

    static void eval(fmpq* q, fmpq const* r, unsigned long int l)
    {
        fmpz_set(fmpq_numref(q), fmpq_numref(r));
        fmpz_mul_ui(fmpq_denref(q), fmpq_denref(r), l);
    }

    static void eval(fmpq* q, unsigned long int l, fmpq const* r)
    {
        if (__FLINTXX_CONSTANT_TRUE(l == 0))
        {
            fmpq_set_ui(q, 0, 1);
        }
        else if (__FLINTXX_CONSTANT_TRUE(l == 1))
        {
            fmpq_inv(q, r);
        }
        else
        {
            eval(q, r, l);
            fmpq_inv(q, q);
        }
    }

    static void eval(fmpq* q, fmpq const* r, signed long int l)
    {
        if (__FLINTXX_CONSTANT_TRUE(l >= 0))
        {
            eval(q, r, static_cast<unsigned long>(l));
        }
        else if (__FLINTXX_CONSTANT_TRUE(l <= 0))
        {
            eval(q, r, -static_cast<unsigned long>(l));
            fmpq_neg(q, q);
        }
        else
        {
            fmpz_set(fmpq_numref(q), fmpq_numref(r));
            fmpz_mul_si(fmpq_denref(q), fmpq_denref(r), l);
        }
    }

    static void eval(fmpq* q, signed long int l, fmpq const* r)
    {
        if (__FLINTXX_CONSTANT_TRUE(l == 0))
        {
            fmpq_set_ui(q, 0, 1);
        }
        else if (__FLINTXX_CONSTANT_TRUE(l == 1))
        {
            fmpq_inv(q, r);
        }
        else if (__FLINTXX_CONSTANT_TRUE(l == -1))
        {
            fmpq_inv(q, r);
            fmpq_neg(q, q);
        }
        else
        {
            eval(q, r, l);
            fmpq_inv(q, q);
        }
    }

    static void eval(fmpq* q, fmpq const* r, double d)
    {
        __flint_binary_multiplies::eval(q, r, static_cast<slong>(1.0 / d));
    }

    static void eval(fmpq* q, double d, fmpq const* r)
    {
        __flint_binary_multiplies::eval(q, static_cast<slong>(1.0 / d), r);
    }
};

struct __flint_binary_modulus
{
    static void eval(fmpz* z, fmpz const* w, fmpz const* v)
    {
        fmpz_fdiv_r(z, w, v);
    }

    static void eval(fmpz* z, fmpz const* w, unsigned long int l)
    {
        ulong tmp = fmpz_tdiv_ui(w, l);
        fmpz_set_ui(z, tmp);
    }

    static void eval(fmpz* z, unsigned long int l, fmpz const* w)
    {
        if (fmpz_sgn(w) >= 0)
        {
            if (fmpz_abs_fits_ui(w))
                fmpz_set_ui(z, l % fmpz_get_ui(w));
            else
                fmpz_set_ui(z, l);
        }
        else
        {
            fmpz_neg(z, w);

            if (fmpz_abs_fits_ui(z))
                fmpz_set_ui(z, l % fmpz_get_ui(z));
            else
                fmpz_set_ui(z, l);
        }
    }

    static void eval(fmpz* z, fmpz const* w, signed long int l)
    {
        fmpz_tdiv_q_ui(z, w, static_cast<ulong>(FLINT_ABS(l)));
    }

    static void eval(fmpz* z, signed long int l, fmpz const* w)
    {
        if (fmpz_abs_fits_ui(w))
        {
            fmpz_set_si(z, l % fmpz_get_si(w));
        }
        else
        {
            fmpz_t tmp;
            fmpz_set_si(tmp, l);
            /* if w is bigger than a long then the remainder is l unchanged,
               unless l==LONG_MIN and w==-LONG_MIN in which case it's 0 */
            fmpz_set_si(z, fmpz_cmpabs(w, tmp) == 0 ? 0 : l);
            fmpz_clear(tmp);
        }
    }

    static void eval(fmpz* z, fmpz const* w, double d)
    {
        eval(z, w, static_cast<slong>(d));
    }

    static void eval(fmpz* z, double d, fmpz const* w)
    {
        eval(z, static_cast<slong>(d), w);
    }
};

struct __flint_binary_and
{
    static void eval(fmpz* z, fmpz const* w, fmpz const* v)
    {
        fmpz_and(z, w, v);
    }

    static void eval(fmpz* z, fmpz const* w, unsigned long int l)
    {
        fmpz_t tmp;
        fmpz_set_ui(tmp, l);
        fmpz_and(z, w, tmp);
        fmpz_clear(tmp);
    }

    static void eval(fmpz* z, unsigned long int l, fmpz const* w)
    {
        eval(z, w, l);
    }

    static void eval(fmpz* z, fmpz const* w, signed long int l)
    {
        fmpz_t tmp;
        fmpz_set_si(tmp, l);
        fmpz_and(z, w, tmp);
        fmpz_clear(tmp);
    }

    static void eval(fmpz* z, signed long int l, fmpz const* w)
    {
        eval(z, w, l);
    }

    static void eval(fmpz* z, fmpz const* w, double d)
    {
        eval(z, w, static_cast<slong>(d));
    }

    static void eval(fmpz* z, double d, fmpz const* w) { eval(z, w, d); }
};

struct __flint_binary_ior
{
    static void eval(fmpz* z, fmpz const* w, fmpz const* v)
    {
        fmpz_or(z, w, v);
    }

    static void eval(fmpz* z, fmpz const* w, unsigned long int l)
    {
        fmpz_t tmp;
        fmpz_set_ui(tmp, l);
        fmpz_or(z, w, tmp);
        fmpz_clear(tmp);
    }

    static void eval(fmpz* z, unsigned long int l, fmpz const* w)
    {
        eval(z, w, l);
    }

    static void eval(fmpz* z, fmpz const* w, signed long int l)
    {
        fmpz_t tmp;
        fmpz_set_si(tmp, l);
        fmpz_or(z, w, tmp);
        fmpz_clear(tmp);
    }

    static void eval(fmpz* z, signed long int l, fmpz const* w)
    {
        eval(z, w, l);
    }

    static void eval(fmpz* z, fmpz const* w, double d)
    {
        eval(z, w, static_cast<slong>(d));
    }

    static void eval(fmpz* z, double d, fmpz const* w) { eval(z, w, d); }
};

struct __flint_binary_xor
{
    static void eval(fmpz* z, fmpz const* w, fmpz const* v)
    {
        fmpz_xor(z, w, v);
    }

    static void eval(fmpz* z, fmpz const* w, unsigned long int l)
    {
        fmpz_t tmp;
        fmpz_set_ui(tmp, l);
        fmpz_xor(z, w, tmp);
        fmpz_clear(tmp);
    }

    static void eval(fmpz* z, unsigned long int l, fmpz const* w)
    {
        eval(z, w, l);
    }

    static void eval(fmpz* z, fmpz const* w, signed long int l)
    {
        fmpz_t tmp;
        fmpz_set_si(tmp, l);
        fmpz_xor(z, w, tmp);
        fmpz_clear(tmp);
    }

    static void eval(fmpz* z, signed long int l, fmpz const* w)
    {
        eval(z, w, l);
    }

    static void eval(fmpz* z, fmpz const* w, double d)
    {
        eval(z, w, static_cast<slong>(d));
    }

    static void eval(fmpz* z, double d, fmpz const* w) { eval(z, w, d); }
};

struct __flint_cmp_function
{
    static int eval(fmpz const* z, fmpz const* w) { return fmpz_cmp(z, w); }

    static int eval(fmpz const* z, unsigned long int l)
    {
        return fmpz_cmp_ui(z, l);
    }

    static int eval(unsigned long int l, fmpz const* z)
    {
        return -fmpz_cmp_ui(z, l);
    }

    static int eval(fmpz const* z, signed long int l)
    {
        return fmpz_cmp_si(z, l);
    }

    static int eval(signed long int l, fmpz const* z)
    {
        return -fmpz_cmp_si(z, l);
    }

    static int eval(fmpz const* z, double d)
    {
        double tmp = fmpz_get_d(z);
        return static_cast<int>(tmp - d);
    }

    static int eval(double d, fmpz const* z) { return -eval(z, d); }

    static int eval(fmpq const* q, fmpq const* r) { return fmpq_cmp(q, r); }

    static int eval(fmpq const* q, unsigned long int l)
    {
        return fmpq_cmp_ui(q, l);
    }

    static int eval(unsigned long int l, fmpq const* q)
    {
        return -fmpq_cmp_ui(q, l);
    }

    static int eval(fmpq const* q, signed long int l)
    {
        return fmpq_cmp_si(q, l);
    }

    static int eval(signed long int l, fmpq const* q)
    {
        return -fmpq_cmp_si(q, l);
    }

    static int eval(fmpq const* q, double d)
    {
        double tmp = fmpq_get_d(q);
        return static_cast<int>(tmp - d);
    }

    static int eval(double d, fmpq const* q) { return -eval(q, d); }

    static int eval(fmpq const* q, fmpz const* z)
    {
        return fmpq_cmp_fmpz(q, z);
    }

    static int eval(fmpz const* z, fmpq const* q)
    {
        return -fmpq_cmp_fmpz(q, z);
    }
};

struct __flint_binary_equal
{
    static bool eval(fmpz const* z, fmpz const* w)
    {
        return fmpz_cmp(z, w) == 0;
    }

    static bool eval(fmpz const* z, unsigned long int l)
    {
        return fmpz_cmp_ui(z, l) == 0;
    }

    static bool eval(unsigned long int l, fmpz const* z) { return eval(z, l); }

    static bool eval(fmpz const* z, signed long int l)
    {
        return fmpz_cmp_si(z, l) == 0;
    }

    static bool eval(signed long int l, fmpz const* z) { return eval(z, l); }

    static bool eval(fmpz const* z, double d)
    {
        return fmpz_cmp_si(z, static_cast<slong>(d)) == 0;
    }

    static bool eval(double d, fmpz const* z) { return eval(z, d); }

    static bool eval(fmpq const* q, fmpq const* r)
    {
        return fmpq_equal(q, r) != 0;
    }

    static bool eval(fmpq const* q, unsigned long int l)
    {
        return ((__FLINTXX_CONSTANT(l) && l == 0)
                || fmpz_cmp_ui(fmpq_denref(q), 1) == 0)
            && fmpz_cmp_ui(fmpq_numref(q), l) == 0;
    }

    static bool eval(unsigned long int l, fmpq const* q) { return eval(q, l); }

    static bool eval(fmpq const* q, signed long int l)
    {
        return ((__FLINTXX_CONSTANT(l) && l == 0)
                || fmpz_cmp_ui(fmpq_denref(q), 1) == 0)
            && fmpz_cmp_si(fmpq_numref(q), l) == 0;
    }

    static bool eval(signed long int l, fmpq const* q) { return eval(q, l); }

    static bool eval(fmpq const* q, double d)
    {
        double tmp = fmpq_get_d(q);
        return tmp == d;
    }

    static bool eval(double d, fmpq const* q) { return eval(q, d); }

    static bool eval(fmpq const* q, fmpz const* z)
    {
        return fmpz_cmp_ui(fmpq_denref(q), 1) == 0
            && fmpz_cmp(fmpq_numref(q), z) == 0;
    }

    static bool eval(fmpz const* z, fmpq const* q) { return eval(q, z); }
};

struct __flint_binary_less
{
    static bool eval(fmpz const* z, fmpz const* w)
    {
        return fmpz_cmp(z, w) < 0;
    }

    static bool eval(fmpz const* z, unsigned long int l)
    {
        return fmpz_cmp_ui(z, l) < 0;
    }

    static bool eval(unsigned long int l, fmpz const* z)
    {
        return fmpz_cmp_ui(z, l) > 0;
    }

    static bool eval(fmpz const* z, signed long int l)
    {
        return fmpz_cmp_si(z, l) < 0;
    }

    static bool eval(signed long int l, fmpz const* z)
    {
        return fmpz_cmp_si(z, l) > 0;
    }

    static bool eval(fmpz const* z, double d)
    {
        double tmp = fmpz_get_d(z);
        return tmp < d;
    }

    static bool eval(double d, fmpz const* z)
    {
        double tmp = fmpz_get_d(z);
        return d < tmp;
    }

    static bool eval(fmpq const* q, fmpq const* r)
    {
        return fmpq_cmp(q, r) < 0;
    }

    static bool eval(fmpq const* q, unsigned long int l)
    {
        return fmpq_cmp_ui(q, l) < 0;
    }

    static bool eval(unsigned long int l, fmpq const* q)
    {
        return fmpq_cmp_ui(q, l) > 0;
    }

    static bool eval(fmpq const* q, signed long int l)
    {
        return fmpq_cmp_si(q, l) < 0;
    }

    static bool eval(signed long int l, fmpq const* q)
    {
        return fmpq_cmp_si(q, l) > 0;
    }

    static bool eval(fmpq const* q, double d)
    {
        double tmp = fmpq_get_d(q);
        return tmp < d;
    }

    static bool eval(double d, fmpq const* q) { return -eval(q, d); }

    static bool eval(fmpq const* q, fmpz const* z)
    {
        return fmpq_cmp_fmpz(q, z) < 0;
    }

    static bool eval(fmpz const* z, fmpq const* q)
    {
        return fmpq_cmp_fmpz(q, z) > 0;
    }
};

struct __flint_binary_greater
{
    template <class T, class U>
    static bool eval(T t, U u)
    {
        return __flint_binary_less::eval(u, t);
    }
};

struct __flint_unary_increment
{
    static void eval(fmpz* z) { fmpz_add_ui(z, z, 1); }

    static void eval(fmpq* q)
    {
        fmpz_add(fmpq_numref(q), fmpq_numref(q), fmpq_denref(q));
    }
};

struct __flint_unary_decrement
{
    static void eval(fmpz* z) { fmpz_sub_ui(z, z, 1); }

    static void eval(fmpq* q)
    {
        fmpz_sub(fmpq_numref(q), fmpq_numref(q), fmpq_denref(q));
    }
};

struct __flint_abs_function
{
    static void eval(fmpz* z, fmpz const* w) { fmpz_abs(z, w); }

    static void eval(fmpq* q, fmpq const* r) { fmpq_abs(q, r); }
};

struct __flint_sqrt_function
{
    static void eval(fmpz* z, fmpz const* w) { fmpz_sqrt(z, w); }
};

struct __flint_sgn_function
{
    static int eval(fmpz const* z) { return fmpz_sgn(z); }

    static int eval(fmpq const* q) { return fmpq_sgn(q); }
};

struct __flint_gcd_function
{
    static void eval(fmpz* z, fmpz const* w, fmpz const* v)
    {
        fmpz_gcd(z, w, v);
    }

    static void eval(fmpz* z, fmpz const* w, unsigned long int l)
    {
        fmpz_gcd_ui(z, w, l);
    }

    static void eval(fmpz* z, unsigned long int l, fmpz const* w)
    {
        eval(z, w, l);
    }

    static void eval(fmpz* z, fmpz const* w, signed long int l)
    {
        fmpz_gcd_ui(z, w, static_cast<ulong>(FLINT_ABS(l)));
    }

    static void eval(fmpz* z, signed long int l, fmpz const* w)
    {
        eval(z, w, l);
    }

    static void eval(fmpz* z, fmpz const* w, double d)
    {
        fmpz_gcd_ui(z, w, static_cast<ulong>(d));
    }

    static void eval(fmpz* z, double d, fmpz const* w) { eval(z, w, d); }
};

struct __flint_lcm_function
{
    static void eval(fmpz* z, fmpz const* w, fmpz const* v)
    {
        fmpz_lcm(z, w, v);
    }

    static void eval(fmpz* z, fmpz const* w, unsigned long int l)
    {
        fmpz_t tmp;
        fmpz_set_ui(tmp, l);
        fmpz_lcm(z, w, tmp);
        fmpz_clear(tmp);
    }

    static void eval(fmpz* z, unsigned long int l, fmpz const* w)
    {
        eval(z, w, l);
    }

    static void eval(fmpz* z, fmpz const* w, signed long int l)
    {
        eval(z, static_cast<ulong>(FLINT_ABS(l)), w);
    }

    static void eval(fmpz* z, signed long int l, fmpz const* w)
    {
        eval(z, w, l);
    }

    static void eval(fmpz* z, fmpz const* w, double d)
    {
        eval(z, w, static_cast<slong>(d));
    }

    static void eval(fmpz* z, double d, fmpz const* w) { eval(z, w, d); }
};
/**************** Auxiliary classes ****************/

// general expression template class
template <class T, class U>
class __flint_expr;

// templates for resolving expression types
template <class T>
struct __flint_resolve_ref
{
    using ref_type = T;
};

template <class T, class U>
struct __flint_resolve_ref<__flint_expr<T, U>>
{
    using ref_type = __flint_expr<T, U> const&;
};

template <class T, class U = T>
struct __flint_resolve_expr;

template <>
struct __flint_resolve_expr<fmpz_t>
{
    using value_type  = fmpz_t;
    using ptr_type    = fmpz*;
    using srcptr_type = fmpz const*;
};

template <>
struct __flint_resolve_expr<fmpq_t>
{
    using value_type  = fmpq_t;
    using ptr_type    = fmpq*;
    using srcptr_type = fmpq const*;
};

template <>
struct __flint_resolve_expr<fmpz_t, fmpq_t>
{
    using value_type = fmpq_t;
};

template <>
struct __flint_resolve_expr<fmpq_t, fmpz_t>
{
    using value_type = fmpq_t;
};

namespace std
{

template <class T, class U, class V, class W>
struct common_type<__flint_expr<T, U>, __flint_expr<V, W>>
{
private:
    using X = typename __flint_resolve_expr<T, V>::value_type;

public:
    using type = __flint_expr<X, X>;
};

template <class T, class U>
struct common_type<__flint_expr<T, U>>
{
    using type = __flint_expr<T, T>;
};

#define __FLINTXX_DECLARE_COMMON_TYPE(typ)      \
    template <class T, class U>                 \
    struct common_type<__flint_expr<T, U>, typ> \
    {                                           \
        typedef __flint_expr<T, T> type;        \
    };                                          \
                                                \
    template <class T, class U>                 \
    struct common_type<typ, __flint_expr<T, U>> \
    {                                           \
        typedef __flint_expr<T, T> type;        \
    }

__FLINTXX_DECLARE_COMMON_TYPE(signed char);
__FLINTXX_DECLARE_COMMON_TYPE(unsigned char);
__FLINTXX_DECLARE_COMMON_TYPE(signed int);
__FLINTXX_DECLARE_COMMON_TYPE(unsigned int);
__FLINTXX_DECLARE_COMMON_TYPE(signed short int);
__FLINTXX_DECLARE_COMMON_TYPE(unsigned short int);
__FLINTXX_DECLARE_COMMON_TYPE(signed long int);
__FLINTXX_DECLARE_COMMON_TYPE(unsigned long int);
__FLINTXX_DECLARE_COMMON_TYPE(float);
__FLINTXX_DECLARE_COMMON_TYPE(double);
#undef __FLINTXX_DECLARE_COMMON_TYPE

}  // namespace std

// classes for evaluating unary and binary expressions
template <class T, class Op>
struct __flint_unary_expr
{
    typename __flint_resolve_ref<T>::ref_type val;

    __flint_unary_expr(T const& v) : val(v) {}

    __flint_unary_expr() = delete;
};

template <class T, class U, class Op>
struct __flint_binary_expr
{
    typename __flint_resolve_ref<T>::ref_type val1;
    typename __flint_resolve_ref<U>::ref_type val2;

    __flint_binary_expr(T const& v1, U const& v2) : val1(v1), val2(v2) {}

    __flint_binary_expr() = delete;
};

/**************** Macros for in-class declarations ****************/
/* This is just repetitive code that is easier to maintain if it's written
   only once */

#define __FLINTXXP_DECLARE_COMPOUND_OPERATOR(fun) \
    template <class T, class U>                   \
    __flint_expr<value_type, value_type>& fun(const __flint_expr<T, U>&);

#define __FLINTXXN_DECLARE_COMPOUND_OPERATOR(fun) \
    __flint_expr& fun(signed char);               \
    __flint_expr& fun(unsigned char);             \
    __flint_expr& fun(signed int);                \
    __flint_expr& fun(unsigned int);              \
    __flint_expr& fun(signed short int);          \
    __flint_expr& fun(unsigned short int);        \
    __flint_expr& fun(signed long int);           \
    __flint_expr& fun(unsigned long int);         \
    __flint_expr& fun(float);                     \
    __flint_expr& fun(double);

#define __FLINTXX_DECLARE_COMPOUND_OPERATOR(fun) \
    __FLINTXXP_DECLARE_COMPOUND_OPERATOR(fun)    \
    __FLINTXXN_DECLARE_COMPOUND_OPERATOR(fun)

#define __FLINTXX_DECLARE_COMPOUND_OPERATOR_UI(fun) \
    __flint_expr& fun(flint_bitcnt_t);

#define __FLINTXX_DECLARE_INCREMENT_OPERATOR(fun) \
    inline __flint_expr& fun();                   \
    inline __flint_expr fun(int);

// clang-format off
#define __FLINTXX_DEFINE_ARITHMETIC_CONSTRUCTORS                                 \
  __flint_expr(signed char c) { init_si(c); }                           \
  __flint_expr(unsigned char c) { init_ui(c); }                         \
  __flint_expr(signed int i) { init_si(i); }                            \
  __flint_expr(unsigned int i) { init_ui(i); }                          \
  __flint_expr(signed short int s) { init_si(s); }                      \
  __flint_expr(unsigned short int s) { init_ui(s); }                    \
  __flint_expr(signed long int l) { init_si(l); }                       \
  __flint_expr(unsigned long int l) { init_ui(l); }                     \
  __flint_expr(float f) { init_d(static_cast<double>(f)); }             \
  __flint_expr(double d) { init_d(d); }
// clang-format on

#define __FLINTXX_DEFINE_ARITHMETIC_ASSIGNMENTS   \
    __flint_expr& operator=(signed char c)        \
    {                                             \
        assign_si(c);                             \
        return *this;                             \
    }                                             \
    __flint_expr& operator=(unsigned char c)      \
    {                                             \
        assign_ui(c);                             \
        return *this;                             \
    }                                             \
    __flint_expr& operator=(signed int i)         \
    {                                             \
        assign_si(i);                             \
        return *this;                             \
    }                                             \
    __flint_expr& operator=(unsigned int i)       \
    {                                             \
        assign_ui(i);                             \
        return *this;                             \
    }                                             \
    __flint_expr& operator=(signed short int s)   \
    {                                             \
        assign_si(s);                             \
        return *this;                             \
    }                                             \
    __flint_expr& operator=(unsigned short int s) \
    {                                             \
        assign_ui(s);                             \
        return *this;                             \
    }                                             \
    __flint_expr& operator=(signed long int l)    \
    {                                             \
        assign_si(l);                             \
        return *this;                             \
    }                                             \
    __flint_expr& operator=(unsigned long int l)  \
    {                                             \
        assign_ui(l);                             \
        return *this;                             \
    }                                             \
    __flint_expr& operator=(float f)              \
    {                                             \
        assign_d(static_cast<double>(f));         \
        return *this;                             \
    }                                             \
    __flint_expr& operator=(double d)             \
    {                                             \
        assign_d(d);                              \
        return *this;                             \
    }

#define __FLINTP_DECLARE_UNARY_STATIC_MEMFUN(T, fun, eval_fun)               \
    template <class U>                                                       \
    static __flint_expr<T, __flint_unary_expr<__flint_expr<T, U>, eval_fun>> \
    fun(const __flint_expr<T, U>& expr);

#define __FLINTNN_DECLARE_UNARY_STATIC_MEMFUN(T, fun, eval_fun, type, bigtype) \
    static inline __flint_expr<T, __flint_unary_expr<(bigtype), eval_fun>>     \
    fun(type expr);

#define __FLINTNS_DECLARE_UNARY_STATIC_MEMFUN(T, fun, eval_fun, type) \
    __FLINTNN_DECLARE_UNARY_STATIC_MEMFUN(T, fun, eval_fun, type, signed long)
#define __FLINTNU_DECLARE_UNARY_STATIC_MEMFUN(T, fun, eval_fun, type) \
    __FLINTNN_DECLARE_UNARY_STATIC_MEMFUN(T, fun, eval_fun, type, unsigned long)
#define __FLINTND_DECLARE_UNARY_STATIC_MEMFUN(T, fun, eval_fun, type) \
    __FLINTNN_DECLARE_UNARY_STATIC_MEMFUN(T, fun, eval_fun, type, double)

#define __FLINTN_DECLARE_UNARY_STATIC_MEMFUN(T, fun, eval_fun)                 \
    __FLINTNS_DECLARE_UNARY_STATIC_MEMFUN(T, fun, eval_fun, signed char)       \
    __FLINTNU_DECLARE_UNARY_STATIC_MEMFUN(T, fun, eval_fun, unsigned char)     \
    __FLINTNS_DECLARE_UNARY_STATIC_MEMFUN(T, fun, eval_fun, signed int)        \
    __FLINTNU_DECLARE_UNARY_STATIC_MEMFUN(T, fun, eval_fun, unsigned int)      \
    __FLINTNS_DECLARE_UNARY_STATIC_MEMFUN(T, fun, eval_fun, signed short int)  \
    __FLINTNU_DECLARE_UNARY_STATIC_MEMFUN(T, fun, eval_fun,                    \
                                          unsigned short int)                  \
    __FLINTNS_DECLARE_UNARY_STATIC_MEMFUN(T, fun, eval_fun, signed long int)   \
    __FLINTNU_DECLARE_UNARY_STATIC_MEMFUN(T, fun, eval_fun, unsigned long int) \
    __FLINTND_DECLARE_UNARY_STATIC_MEMFUN(T, fun, eval_fun, float)             \
    __FLINTND_DECLARE_UNARY_STATIC_MEMFUN(T, fun, eval_fun, double)

#define __FLINT_DECLARE_UNARY_STATIC_MEMFUN(T, fun, eval_fun) \
    __FLINTP_DECLARE_UNARY_STATIC_MEMFUN(T, fun, eval_fun)    \
    __FLINTN_DECLARE_UNARY_STATIC_MEMFUN(T, fun, eval_fun)

/**************** fmpz_class -- wrapper for fmpz_t ****************/

template <>
class __flint_expr<fmpz_t, fmpz_t>
{
    using value_type = fmpz_t;

    fmpz mp{};

    // Helper functions used for all arithmetic types
    void assign_ui(unsigned long l) { fmpz_set_ui(&mp, l); }

    void assign_si(signed long l) { fmpz_set_si(&mp, l); }

    void assign_d(double d) { fmpz_set_d(&mp, d); }

    void init_ui(unsigned long l) { fmpz_init_set_ui(&mp, l); }

    void init_si(signed long l) { fmpz_init_set_si(&mp, l); }

    void init_d(double d)
    {
        fmpz_init(&mp);
        fmpz_set_d(&mp, d);
    }

public:
    // constructors and destructor
    __flint_expr() noexcept { fmpz_init(&mp); }

    __flint_expr(__flint_expr const& z) { fmpz_init_set(&mp, &z.mp); }

    __flint_expr(__flint_expr&& z) noexcept : mp(z.mp) { fmpz_init(&z.mp); }

    template <class T>
    __flint_expr(__flint_expr<fmpz_t, T> const& expr)
    {
        fmpz_init(&mp);
        __flint_set_expr(&mp, expr);
    }

    template <class T, class U>
    explicit __flint_expr(__flint_expr<T, U> const& expr)
    {
        fmpz_init(&mp);
        __flint_set_expr(&mp, expr);
    }

    __FLINTXX_DEFINE_ARITHMETIC_CONSTRUCTORS

    explicit __flint_expr(char const* s, int base = 0)
    {
        fmpz_init(&mp);

        if (fmpz_set_str(&mp, s, base) != 0)
        {
            fmpz_clear(&mp);
            throw std::invalid_argument("fmpz_set_str");
        }
    }

    explicit __flint_expr(std::string const& s, int base = 0)
    {
        fmpz_init(&mp);

        if (fmpz_set_str(&mp, s.c_str(), base) != 0)
        {
            fmpz_clear(&mp);
            throw std::invalid_argument("fmpz_set_str");
        }
    }

    explicit __flint_expr(fmpz const* z) { fmpz_init_set(&mp, z); }

    ~__flint_expr() { fmpz_clear(&mp); }

    void swap(__flint_expr& z) noexcept { std::swap(mp, z.mp); }

    // assignment operators
    __flint_expr& operator=(__flint_expr const& z)
    {
        if (&z == this)
            return *this;

        fmpz_set(&mp, &z.mp);
        return *this;
    }

    __flint_expr& operator=(__flint_expr&& z) noexcept
    {
        swap(z);
        return *this;
    }

    template <class T, class U>
    __flint_expr& operator=(__flint_expr<T, U> const& expr)
    {
        __flint_set_expr(&mp, expr);
        return *this;
    }

    __FLINTXX_DEFINE_ARITHMETIC_ASSIGNMENTS

    __flint_expr& operator=(char const* s)
    {
        if (fmpz_set_str(&mp, s, 0) != 0)
            throw std::invalid_argument("fmpz_set_str");

        return *this;
    }

    __flint_expr& operator=(std::string const& s)
    {
        if (fmpz_set_str(&mp, s.c_str(), 0) != 0)
            throw std::invalid_argument("fmpz_set_str");

        return *this;
    }

    // string input/output functions
    int set_str(char const* s, int base) { return fmpz_set_str(&mp, s, base); }

    int set_str(std::string const& s, int base)
    {
        return fmpz_set_str(&mp, s.c_str(), base);
    }

    std::string get_str(int base = 10) const
    {
        char* str = fmpz_get_str(nullptr, base, &mp);
        std::string tmp{str};
        flint_free(str);

        return tmp;
    }

    // conversion functions
    fmpz const* __get_mp() const { return &mp; }  // NOLINT

    fmpz* __get_mp() { return &mp; }  // NOLINT

    fmpz const* get_fmpz_t() const { return &mp; }

    fmpz* get_fmpz_t() { return &mp; }

    signed long int get_si() const { return fmpz_get_si(&mp); }  // NOLINT

    unsigned long int get_ui() const { return fmpz_get_ui(&mp); }  // NOLINT

    double get_d() const { return fmpz_get_d(&mp); }  // NOLINT

    explicit operator bool() const { return not fmpz_is_zero(&mp); }

    // member operators
    __FLINTXX_DECLARE_COMPOUND_OPERATOR(operator+=)
    __FLINTXX_DECLARE_COMPOUND_OPERATOR(operator-=)
    __FLINTXX_DECLARE_COMPOUND_OPERATOR(operator*=)
    __FLINTXX_DECLARE_COMPOUND_OPERATOR(operator/=)
    __FLINTXX_DECLARE_COMPOUND_OPERATOR(operator%=)

    __FLINTXX_DECLARE_COMPOUND_OPERATOR(operator&=)
    __FLINTXX_DECLARE_COMPOUND_OPERATOR(operator|=)
    __FLINTXX_DECLARE_COMPOUND_OPERATOR(operator^=)

    __FLINTXX_DECLARE_COMPOUND_OPERATOR_UI(operator<<=)
    __FLINTXX_DECLARE_COMPOUND_OPERATOR_UI(operator>>=)

    __FLINTXX_DECLARE_INCREMENT_OPERATOR(operator++)  // NOLINT
    __FLINTXX_DECLARE_INCREMENT_OPERATOR(operator--)  // NOLINT
};

using fmpz_class = __flint_expr<fmpz_t, fmpz_t>;

/**************** fmpq_class -- wrapper for fmpq_t ****************/

template <>
class __flint_expr<fmpq_t, fmpq_t>
{
    using value_type = fmpq_t;

    fmpq mp{};

    // Helper functions used for all arithmetic types
    void assign_ui(unsigned long l) { fmpq_set_ui(&mp, l, 1); }  // NOLINT

    void assign_si(signed long l) { fmpq_set_si(&mp, l, 1); }  // NOLINT

    void assign_d(double d)
    {
        mpq_t tmp;
        mpq_init(tmp);

        mpq_set_d(tmp, d);

        fmpq_set_mpq(&mp, tmp);

        mpq_clear(tmp);
    }

    void init_ui(unsigned long l)  // NOLINT
    {
        fmpq_init(&mp);
        get_num() = l;
    }

    void init_si(signed long l)  // NOLINT
    {
        fmpq_init(&mp);
        get_num() = l;
    }

    void init_d(double d)
    {
        fmpq_init(&mp);
        assign_d(d);
    }

public:
    void canonicalize() { fmpq_canonicalise(&mp); }

    // constructors and destructor
    __flint_expr() { fmpq_init(&mp); }

    __flint_expr(__flint_expr const& q)
    {
        fmpz_init_set(fmpq_numref(&mp), fmpq_numref(&q.mp));
        fmpz_init_set(fmpq_denref(&mp), fmpq_denref(&q.mp));
    }

    __flint_expr(__flint_expr&& q) noexcept : mp(q.mp) { fmpq_init(&q.mp); }

    template <class T>
    __flint_expr(__flint_expr<fmpz_t, T> const& expr)
    {
        fmpq_init(&mp);
        __flint_set_expr(&mp, expr);
    }

    template <class T>
    __flint_expr(__flint_expr<fmpq_t, T> const& expr)
    {
        fmpq_init(&mp);
        __flint_set_expr(&mp, expr);
    }

    template <class T, class U>
    explicit __flint_expr(__flint_expr<T, U> const& expr)
    {
        fmpq_init(&mp);
        __flint_set_expr(&mp, expr);
    }

    // conversion functions

    // casting a reference to an fmpz to fmpz_class & is a dirty hack,
    // but works because the internal representation of fmpz_class is
    // exactly an fmpz
    fmpz_class const& get_num() const
    {
        return reinterpret_cast<fmpz_class const&>(*fmpq_numref(&mp));
    }

    fmpz_class& get_num()
    {
        return reinterpret_cast<fmpz_class&>(*fmpq_numref(&mp));
    }

    fmpz_class const& get_den() const
    {
        return reinterpret_cast<fmpz_class const&>(*fmpq_denref(&mp));
    }

    fmpz_class& get_den()
    {
        return reinterpret_cast<fmpz_class&>(*fmpq_denref(&mp));
    }

    __FLINTXX_DEFINE_ARITHMETIC_CONSTRUCTORS

    explicit __flint_expr(char const* s, int base = 0)
    {
        fmpq_init(&mp);
        // If s is the literal 0, we meant to call another constructor.
        // If s just happens to evaluate to 0, we would crash, so whatever.
        if (s == nullptr)
        {
            // Don't turn fmpq_class(0,0) into 0
            fmpz_set_si(fmpq_denref(&mp), base);
        }
        else if (fmpq_set_str(&mp, s, base) != 0)
        {
            fmpq_clear(&mp);

            throw std::invalid_argument("fmpq_set_str");
        }
    }

    explicit __flint_expr(std::string const& s, int base = 0)
    {
        fmpq_init(&mp);

        if (fmpq_set_str(&mp, s.c_str(), base) != 0)
        {
            fmpq_clear(&mp);

            throw std::invalid_argument("fmpq_set_str");
        }
    }

    explicit __flint_expr(fmpq const* q)
    {
        fmpz_init_set(fmpq_numref(&mp), fmpq_numref(q));
        fmpz_init_set(fmpq_denref(&mp), fmpq_denref(q));
    }

    __flint_expr(fmpz_class const& num, fmpz_class const& den)
    {
        fmpz_init_set(fmpq_numref(&mp), num.get_fmpz_t());
        fmpz_init_set(fmpq_denref(&mp), den.get_fmpz_t());
    }

    ~__flint_expr() { fmpq_clear(&mp); }

    void swap(__flint_expr& q) noexcept { std::swap(mp, q.mp); }

    // assignment operators
    __flint_expr& operator=(__flint_expr const& q)
    {
        if (&q == this)
            return *this;

        fmpq_set(&mp, &q.mp);
        return *this;
    }

    __flint_expr& operator=(__flint_expr&& q) noexcept
    {
        swap(q);
        return *this;
    }

    __flint_expr& operator=(fmpz_class&& z) noexcept
    {
        get_num() = std::move(z);
        get_den() = 1U;
        return *this;
    }

    template <class T, class U>
    __flint_expr& operator=(__flint_expr<T, U> const& expr)
    {
        __flint_set_expr(&mp, expr);
        return *this;
    }

    __FLINTXX_DEFINE_ARITHMETIC_ASSIGNMENTS

    __flint_expr& operator=(char const* s)
    {
        if (fmpq_set_str(&mp, s, 0) != 0)
            throw std::invalid_argument("fmpq_set_str");

        return *this;
    }

    __flint_expr& operator=(std::string const& s)
    {
        if (fmpq_set_str(&mp, s.c_str(), 0) != 0)
            throw std::invalid_argument("fmpq_set_str");

        return *this;
    }

    // string input/output functions
    int set_str(char const* s, int base) { return fmpq_set_str(&mp, s, base); }

    int set_str(std::string const& s, int base)
    {
        return fmpq_set_str(&mp, s.c_str(), base);
    }

    std::string get_str(int base = 10) const
    {
        char* str = fmpq_get_str(nullptr, base, &mp);
        std::string tmp{str};
        flint_free(str);

        return tmp;
    }

    fmpq const* __get_mp() const { return &mp; }  // NOLINT

    fmpq* __get_mp() { return &mp; }  // NOLINT

    fmpq const* get_fmpq_t() const { return &mp; }

    fmpq* get_fmpq_t() { return &mp; }

    fmpz const* get_num_fmpz_t() const { return fmpq_numref(&mp); }

    fmpz* get_num_fmpz_t() { return fmpq_numref(&mp); }

    fmpz const* get_den_fmpz_t() const { return fmpq_denref(&mp); }

    fmpz* get_den_fmpz_t() { return fmpq_denref(&mp); }

    double get_d() const { return fmpq_get_d(&mp); }

    explicit operator bool() const { return not fmpq_is_zero(&mp); }

    // compound assignments
    __FLINTXX_DECLARE_COMPOUND_OPERATOR(operator+=)
    __FLINTXX_DECLARE_COMPOUND_OPERATOR(operator-=)
    __FLINTXX_DECLARE_COMPOUND_OPERATOR(operator*=)
    __FLINTXX_DECLARE_COMPOUND_OPERATOR(operator/=)

    __FLINTXX_DECLARE_COMPOUND_OPERATOR_UI(operator<<=)
    __FLINTXX_DECLARE_COMPOUND_OPERATOR_UI(operator>>=)

    __FLINTXX_DECLARE_INCREMENT_OPERATOR(operator++)  // NOLINT
    __FLINTXX_DECLARE_INCREMENT_OPERATOR(operator--)  // NOLINT
};

using fmpq_class = __flint_expr<fmpq_t, fmpq_t>;

/**************** User-defined literals ****************/

inline fmpz_class operator""_fmpz(char const* s)
{
    return fmpz_class(s);
}

inline fmpq_class operator""_fmpq(char const* s)
{
    fmpq_class q;
    q.get_num() = s;
    return q;
}

/**************** Functions for type conversion ****************/

inline void __flint_set_expr(fmpz* z, fmpz_class const& w)
{
    fmpz_set(z, w.get_fmpz_t());
}

template <class T>
inline void __flint_set_expr(fmpz* z, __flint_expr<fmpz_t, T> const& expr)
{
    expr.eval(z);
}

template <class T>
inline void __flint_set_expr(fmpz* z, __flint_expr<fmpq_t, T> const& expr)
{
    fmpq_class const& temp(expr);
    fmpz_set(z, temp.get_num_fmpz_t());
}

inline void __flint_set_expr(fmpq* q, fmpz_class const& z)
{
    fmpq_set_fmpz(q, z.get_fmpz_t());
}

template <class T>
inline void __flint_set_expr(fmpq* q, __flint_expr<fmpz_t, T> const& expr)
{
    __flint_set_expr(fmpq_numref(q), expr);
    fmpz_set_ui(fmpq_denref(q), 1);
}

inline void __flint_set_expr(fmpq* q, fmpq_class const& r)
{
    fmpq_set(q, r.get_fmpq_t());
}

template <class T>
inline void __flint_set_expr(fmpq* q, __flint_expr<fmpq_t, T> const& expr)
{
    expr.eval(q);
}

/* Temporary objects */

template <class T>
class __flint_temp
{
    __flint_expr<T, T> val;

public:
    template <class U, class V>
    __flint_temp(U const& u, V /*unused*/) : val(u)
    {}

    typename __flint_resolve_expr<T>::srcptr_type __get_mp() const
    {
        return val.__get_mp();
    }
};

/**************** Specializations of __flint_expr ****************/
/* The eval() method of __flint_expr<T, U> evaluates the corresponding
   expression and assigns the result to its argument, which is either an
   fmpz_t or fmpq_t as specified by the T argument.
   Compound expressions are evaluated recursively (temporaries are created
   to hold intermediate values), while for simple expressions the eval()
   method of the appropriate function object (available as the Op argument
   of either __flint_unary_expr<T, Op> or __flint_binary_expr<T, U, Op>) is
   called. */

/**************** Unary expressions ****************/
/* cases:
   - simple:   argument is mp*_class, that is, __flint_expr<T, T>
   - compound: argument is __flint_expr<T, U> (with U not equal to T) */

// simple expressions

template <class T, class Op>
class __flint_expr<T, __flint_unary_expr<__flint_expr<T, T>, Op>>
{
    using val_type = __flint_expr<T, T>;

    __flint_unary_expr<val_type, Op> expr;

public:
    explicit __flint_expr(val_type const& val) : expr(val) {}

    void eval(typename __flint_resolve_expr<T>::ptr_type p) const
    {
        Op::eval(p, expr.val.__get_mp());
    }

    val_type const& get_val() const { return expr.val; }
};

// simple expressions, U is a built-in numerical type

template <class T, class U, class Op>
class __flint_expr<T, __flint_unary_expr<U, Op>>
{
    using val_type = U;

    __flint_unary_expr<val_type, Op> expr;

public:
    explicit __flint_expr(val_type const& val) : expr(val) {}

    void eval(typename __flint_resolve_expr<T>::ptr_type p) const
    {
        Op::eval(p, expr.val);
    }

    val_type const& get_val() const { return expr.val; }
};

// compound expressions

template <class T, class U, class Op>
class __flint_expr<T, __flint_unary_expr<__flint_expr<T, U>, Op>>
{
    using val_type = __flint_expr<T, U>;

    __flint_unary_expr<val_type, Op> expr;

public:
    explicit __flint_expr(val_type const& val) : expr(val) {}

    void eval(typename __flint_resolve_expr<T>::ptr_type p) const
    {
        expr.val.eval(p);
        Op::eval(p, p);
    }

    val_type const& get_val() const { return expr.val; }
};

/**************** Binary expressions ****************/
/* simple:
   - arguments are both mp*_class
   - one argument is mp*_class, one is a built-in type
   compound:
   - one is mp*_class, one is __flint_expr<T, U>
   - one is __flint_expr<T, U>, one is built-in
   - both arguments are __flint_expr<...> */

// simple expressions

template <class T, class Op>
class __flint_expr<
    T,
    __flint_binary_expr<__flint_expr<T, T>, __flint_expr<T, T>, Op>>
{
    using val1_type = __flint_expr<T, T>;
    using val2_type = __flint_expr<T, T>;

    __flint_binary_expr<val1_type, val2_type, Op> expr;

public:
    __flint_expr(val1_type const& val1, val2_type const& val2) :
            expr(val1, val2)
    {}

    void eval(typename __flint_resolve_expr<T>::ptr_type p) const
    {
        Op::eval(p, expr.val1.__get_mp(), expr.val2.__get_mp());
    }

    val1_type const& get_val1() const { return expr.val1; }

    val2_type const& get_val2() const { return expr.val2; }
};

// simple expressions, U is a built-in numerical type

template <class T, class U, class Op>
class __flint_expr<T, __flint_binary_expr<__flint_expr<T, T>, U, Op>>
{
    using val1_type = __flint_expr<T, T>;
    using val2_type = U;

    __flint_binary_expr<val1_type, val2_type, Op> expr;

public:
    __flint_expr(val1_type const& val1, val2_type const& val2) :
            expr(val1, val2)
    {}

    void eval(typename __flint_resolve_expr<T>::ptr_type p) const
    {
        Op::eval(p, expr.val1.__get_mp(), expr.val2);
    }

    val1_type const& get_val1() const { return expr.val1; }

    val2_type const& get_val2() const { return expr.val2; }
};

template <class T, class U, class Op>
class __flint_expr<T, __flint_binary_expr<U, __flint_expr<T, T>, Op>>
{
    using val1_type = U;
    using val2_type = __flint_expr<T, T>;

    __flint_binary_expr<val1_type, val2_type, Op> expr;

public:
    __flint_expr(val1_type const& val1, val2_type const& val2) :
            expr(val1, val2)
    {}

    void eval(typename __flint_resolve_expr<T>::ptr_type p) const
    {
        Op::eval(p, expr.val1, expr.val2.__get_mp());
    }

    val1_type const& get_val1() const { return expr.val1; }

    val2_type const& get_val2() const { return expr.val2; }
};

// compound expressions, one argument is a subexpression

template <class T, class U, class V, class Op>
class __flint_expr<
    T,
    __flint_binary_expr<__flint_expr<T, T>, __flint_expr<U, V>, Op>>
{
    using val1_type = __flint_expr<T, T>;
    using val2_type = __flint_expr<U, V>;

    __flint_binary_expr<val1_type, val2_type, Op> expr;

public:
    __flint_expr(val1_type const& val1, val2_type const& val2) :
            expr(val1, val2)
    {}

    void eval(typename __flint_resolve_expr<T>::ptr_type p) const
    {
        if (p != expr.val1.__get_mp())
        {
            __flint_set_expr(p, expr.val2);
            Op::eval(p, expr.val1.__get_mp(), p);
        }
        else
        {
            __flint_temp<T> temp(expr.val2, p);
            Op::eval(p, expr.val1.__get_mp(), temp.__get_mp());
        }
    }

    val1_type const& get_val1() const { return expr.val1; }

    val2_type const& get_val2() const { return expr.val2; }
};

template <class T, class U, class V, class Op>
class __flint_expr<
    T,
    __flint_binary_expr<__flint_expr<U, V>, __flint_expr<T, T>, Op>>
{
    using val1_type = __flint_expr<U, V>;
    using val2_type = __flint_expr<T, T>;

    __flint_binary_expr<val1_type, val2_type, Op> expr;

public:
    __flint_expr(val1_type const& val1, val2_type const& val2) :
            expr(val1, val2)
    {}

    void eval(typename __flint_resolve_expr<T>::ptr_type p) const
    {
        if (p != expr.val2.__get_mp())
        {
            __flint_set_expr(p, expr.val1);
            Op::eval(p, p, expr.val2.__get_mp());
        }
        else
        {
            __flint_temp<T> temp(expr.val1, p);
            Op::eval(p, temp.__get_mp(), expr.val2.__get_mp());
        }
    }

    val1_type const& get_val1() const { return expr.val1; }

    val2_type const& get_val2() const { return expr.val2; }
};

template <class T, class U, class Op>
class __flint_expr<
    T,
    __flint_binary_expr<__flint_expr<T, T>, __flint_expr<T, U>, Op>>
{
    using val1_type = __flint_expr<T, T>;
    using val2_type = __flint_expr<T, U>;

    __flint_binary_expr<val1_type, val2_type, Op> expr;

public:
    __flint_expr(val1_type const& val1, val2_type const& val2) :
            expr(val1, val2)
    {}

    void eval(typename __flint_resolve_expr<T>::ptr_type p) const
    {
        if (p != expr.val1.__get_mp())
        {
            __flint_set_expr(p, expr.val2);
            Op::eval(p, expr.val1.__get_mp(), p);
        }
        else
        {
            __flint_temp<T> temp(expr.val2, p);
            Op::eval(p, expr.val1.__get_mp(), temp.__get_mp());
        }
    }

    val1_type const& get_val1() const { return expr.val1; }

    val2_type const& get_val2() const { return expr.val2; }
};

template <class T, class U, class Op>
class __flint_expr<
    T,
    __flint_binary_expr<__flint_expr<T, U>, __flint_expr<T, T>, Op>>
{
    using val1_type = __flint_expr<T, U>;
    using val2_type = __flint_expr<T, T>;

    __flint_binary_expr<val1_type, val2_type, Op> expr;

public:
    __flint_expr(val1_type const& val1, val2_type const& val2) :
            expr(val1, val2)
    {}

    void eval(typename __flint_resolve_expr<T>::ptr_type p) const
    {
        if (p != expr.val2.__get_mp())
        {
            __flint_set_expr(p, expr.val1);
            Op::eval(p, p, expr.val2.__get_mp());
        }
        else
        {
            __flint_temp<T> temp(expr.val1, p);
            Op::eval(p, temp.__get_mp(), expr.val2.__get_mp());
        }
    }

    val1_type const& get_val1() const { return expr.val1; }

    val2_type const& get_val2() const { return expr.val2; }
};

// one argument is a subexpression, one is a built-in

template <class T, class U, class V, class Op>
class __flint_expr<T, __flint_binary_expr<__flint_expr<T, U>, V, Op>>
{
    using val1_type = __flint_expr<T, U>;
    using val2_type = V;

    __flint_binary_expr<val1_type, val2_type, Op> expr;

public:
    __flint_expr(val1_type const& val1, val2_type const& val2) :
            expr(val1, val2)
    {}

    void eval(typename __flint_resolve_expr<T>::ptr_type p) const
    {
        expr.val1.eval(p);
        Op::eval(p, p, expr.val2);
    }

    val1_type const& get_val1() const { return expr.val1; }

    val2_type const& get_val2() const { return expr.val2; }
};

template <class T, class U, class V, class Op>
class __flint_expr<T, __flint_binary_expr<U, __flint_expr<T, V>, Op>>
{
    using val1_type = U;
    using val2_type = __flint_expr<T, V>;

    __flint_binary_expr<val1_type, val2_type, Op> expr;

public:
    __flint_expr(val1_type const& val1, val2_type const& val2) :
            expr(val1, val2)
    {}

    void eval(typename __flint_resolve_expr<T>::ptr_type p) const
    {
        expr.val2.eval(p);
        Op::eval(p, expr.val1, p);
    }

    val1_type const& get_val1() const { return expr.val1; }

    val2_type const& get_val2() const { return expr.val2; }
};

// both arguments are subexpressions

template <class T, class U, class V, class W, class Op>
class __flint_expr<
    T,
    __flint_binary_expr<__flint_expr<T, U>, __flint_expr<V, W>, Op>>
{
    using val1_type = __flint_expr<T, U>;
    using val2_type = __flint_expr<V, W>;

    __flint_binary_expr<val1_type, val2_type, Op> expr;

public:
    __flint_expr(val1_type const& val1, val2_type const& val2) :
            expr(val1, val2)
    {}

    void eval(typename __flint_resolve_expr<T>::ptr_type p) const
    {
        __flint_temp<T> temp2(expr.val2, p);
        expr.val1.eval(p);
        Op::eval(p, p, temp2.__get_mp());
    }

    val1_type const& get_val1() const { return expr.val1; }

    val2_type const& get_val2() const { return expr.val2; }
};

template <class T, class U, class V, class W, class Op>
class __flint_expr<
    T,
    __flint_binary_expr<__flint_expr<U, V>, __flint_expr<T, W>, Op>>
{
    using val1_type = __flint_expr<U, V>;
    using val2_type = __flint_expr<T, W>;

    __flint_binary_expr<val1_type, val2_type, Op> expr;

public:
    __flint_expr(val1_type const& val1, val2_type const& val2) :
            expr(val1, val2)
    {}

    void eval(typename __flint_resolve_expr<T>::ptr_type p) const
    {
        __flint_temp<T> temp1(expr.val1, p);
        expr.val2.eval(p);
        Op::eval(p, temp1.__get_mp(), p);
    }

    val1_type const& get_val1() const { return expr.val1; }

    val2_type const& get_val2() const { return expr.val2; }
};

template <class T, class U, class V, class Op>
class __flint_expr<
    T,
    __flint_binary_expr<__flint_expr<T, U>, __flint_expr<T, V>, Op>>
{
    using val1_type = __flint_expr<T, U>;
    using val2_type = __flint_expr<T, V>;

    __flint_binary_expr<val1_type, val2_type, Op> expr;

public:
    __flint_expr(val1_type const& val1, val2_type const& val2) :
            expr(val1, val2)
    {}

    void eval(typename __flint_resolve_expr<T>::ptr_type p) const
    {
        __flint_temp<T> temp2(expr.val2, p);
        expr.val1.eval(p);
        Op::eval(p, p, temp2.__get_mp());
    }

    val1_type const& get_val1() const { return expr.val1; }

    val2_type const& get_val2() const { return expr.val2; }
};

/**************** Special cases ****************/

/* Some operations (i.e., add and subtract) with mixed mpz/mpq arguments
   can be done directly without first converting the mpz to mpq.
   Appropriate specializations of __flint_expr are required. */

#define __FLINTZQ_DEFINE_EXPR(eval_fun)                                        \
                                                                               \
    template <>                                                                \
    class __flint_expr<fmpq_t,                                                 \
                       __flint_binary_expr<fmpz_class, fmpq_class, eval_fun>>  \
    {                                                                          \
    private:                                                                   \
        typedef fmpz_class val1_type;                                          \
        typedef fmpq_class val2_type;                                          \
                                                                               \
        __flint_binary_expr<val1_type, val2_type, eval_fun> expr;              \
                                                                               \
    public:                                                                    \
        __flint_expr(const val1_type& val1, const val2_type& val2) :           \
                expr(val1, val2)                                               \
        {}                                                                     \
        void eval(fmpq* q) const                                               \
        {                                                                      \
            eval_fun::eval(q, expr.val1.get_fmpz_t(), expr.val2.get_fmpq_t()); \
        }                                                                      \
        const val1_type& get_val1() const                                      \
        {                                                                      \
            return expr.val1;                                                  \
        }                                                                      \
        const val2_type& get_val2() const                                      \
        {                                                                      \
            return expr.val2;                                                  \
        }                                                                      \
    };                                                                         \
                                                                               \
    template <>                                                                \
    class __flint_expr<fmpq_t,                                                 \
                       __flint_binary_expr<fmpq_class, fmpz_class, eval_fun>>  \
    {                                                                          \
    private:                                                                   \
        typedef fmpq_class val1_type;                                          \
        typedef fmpz_class val2_type;                                          \
                                                                               \
        __flint_binary_expr<val1_type, val2_type, eval_fun> expr;              \
                                                                               \
    public:                                                                    \
        __flint_expr(const val1_type& val1, const val2_type& val2) :           \
                expr(val1, val2)                                               \
        {}                                                                     \
        void eval(fmpq* q) const                                               \
        {                                                                      \
            eval_fun::eval(q, expr.val1.get_fmpq_t(), expr.val2.get_fmpz_t()); \
        }                                                                      \
        const val1_type& get_val1() const                                      \
        {                                                                      \
            return expr.val1;                                                  \
        }                                                                      \
        const val2_type& get_val2() const                                      \
        {                                                                      \
            return expr.val2;                                                  \
        }                                                                      \
    };                                                                         \
                                                                               \
    template <class T>                                                         \
    class __flint_expr<                                                        \
        fmpq_t,                                                                \
        __flint_binary_expr<fmpz_class, __flint_expr<fmpq_t, T>, eval_fun>>    \
    {                                                                          \
    private:                                                                   \
        typedef fmpz_class val1_type;                                          \
        typedef __flint_expr<fmpq_t, T> val2_type;                             \
                                                                               \
        __flint_binary_expr<val1_type, val2_type, eval_fun> expr;              \
                                                                               \
    public:                                                                    \
        __flint_expr(const val1_type& val1, const val2_type& val2) :           \
                expr(val1, val2)                                               \
        {}                                                                     \
        void eval(fmpq* q) const                                               \
        {                                                                      \
            fmpq_class temp(expr.val2);                                        \
            eval_fun::eval(q, expr.val1.get_fmpz_t(), temp.get_fmpq_t());      \
        }                                                                      \
        const val1_type& get_val1() const                                      \
        {                                                                      \
            return expr.val1;                                                  \
        }                                                                      \
        const val2_type& get_val2() const                                      \
        {                                                                      \
            return expr.val2;                                                  \
        }                                                                      \
    };                                                                         \
                                                                               \
    template <class T>                                                         \
    class __flint_expr<                                                        \
        fmpq_t,                                                                \
        __flint_binary_expr<fmpq_class, __flint_expr<fmpz_t, T>, eval_fun>>    \
    {                                                                          \
    private:                                                                   \
        typedef fmpq_class val1_type;                                          \
        typedef __flint_expr<fmpz_t, T> val2_type;                             \
                                                                               \
        __flint_binary_expr<val1_type, val2_type, eval_fun> expr;              \
                                                                               \
    public:                                                                    \
        __flint_expr(const val1_type& val1, const val2_type& val2) :           \
                expr(val1, val2)                                               \
        {}                                                                     \
        void eval(fmpq* q) const                                               \
        {                                                                      \
            fmpz_class temp(expr.val2);                                        \
            eval_fun::eval(q, expr.val1.get_fmpq_t(), temp.get_fmpz_t());      \
        }                                                                      \
        const val1_type& get_val1() const                                      \
        {                                                                      \
            return expr.val1;                                                  \
        }                                                                      \
        const val2_type& get_val2() const                                      \
        {                                                                      \
            return expr.val2;                                                  \
        }                                                                      \
    };                                                                         \
                                                                               \
    template <class T>                                                         \
    class __flint_expr<fmpq_t, __flint_binary_expr<__flint_expr<fmpz_t, T>,    \
                                                   fmpq_class, eval_fun>>      \
    {                                                                          \
    private:                                                                   \
        typedef __flint_expr<fmpz_t, T> val1_type;                             \
        typedef fmpq_class val2_type;                                          \
                                                                               \
        __flint_binary_expr<val1_type, val2_type, eval_fun> expr;              \
                                                                               \
    public:                                                                    \
        __flint_expr(const val1_type& val1, const val2_type& val2) :           \
                expr(val1, val2)                                               \
        {}                                                                     \
        void eval(fmpq* q) const                                               \
        {                                                                      \
            fmpz_class temp(expr.val1);                                        \
            eval_fun::eval(q, temp.get_fmpz_t(), expr.val2.get_fmpq_t());      \
        }                                                                      \
        const val1_type& get_val1() const                                      \
        {                                                                      \
            return expr.val1;                                                  \
        }                                                                      \
        const val2_type& get_val2() const                                      \
        {                                                                      \
            return expr.val2;                                                  \
        }                                                                      \
    };                                                                         \
                                                                               \
    template <class T>                                                         \
    class __flint_expr<fmpq_t, __flint_binary_expr<__flint_expr<fmpq_t, T>,    \
                                                   fmpz_class, eval_fun>>      \
    {                                                                          \
    private:                                                                   \
        typedef __flint_expr<fmpq_t, T> val1_type;                             \
        typedef fmpz_class val2_type;                                          \
                                                                               \
        __flint_binary_expr<val1_type, val2_type, eval_fun> expr;              \
                                                                               \
    public:                                                                    \
        __flint_expr(const val1_type& val1, const val2_type& val2) :           \
                expr(val1, val2)                                               \
        {}                                                                     \
        void eval(fmpq* q) const                                               \
        {                                                                      \
            fmpq_class temp(expr.val1);                                        \
            eval_fun::eval(q, temp.get_fmpq_t(), expr.val2.get_fmpz_t());      \
        }                                                                      \
        const val1_type& get_val1() const                                      \
        {                                                                      \
            return expr.val1;                                                  \
        }                                                                      \
        const val2_type& get_val2() const                                      \
        {                                                                      \
            return expr.val2;                                                  \
        }                                                                      \
    };                                                                         \
                                                                               \
    template <class T, class U>                                                \
    class __flint_expr<fmpq_t,                                                 \
                       __flint_binary_expr<__flint_expr<fmpz_t, T>,            \
                                           __flint_expr<fmpq_t, U>, eval_fun>> \
    {                                                                          \
    private:                                                                   \
        typedef __flint_expr<fmpz_t, T> val1_type;                             \
        typedef __flint_expr<fmpq_t, U> val2_type;                             \
                                                                               \
        __flint_binary_expr<val1_type, val2_type, eval_fun> expr;              \
                                                                               \
    public:                                                                    \
        __flint_expr(const val1_type& val1, const val2_type& val2) :           \
                expr(val1, val2)                                               \
        {}                                                                     \
        void eval(fmpq* q) const                                               \
        {                                                                      \
            fmpz_class temp1(expr.val1);                                       \
            expr.val2.eval(q);                                                 \
            eval_fun::eval(q, temp1.get_fmpz_t(), q);                          \
        }                                                                      \
        const val1_type& get_val1() const                                      \
        {                                                                      \
            return expr.val1;                                                  \
        }                                                                      \
        const val2_type& get_val2() const                                      \
        {                                                                      \
            return expr.val2;                                                  \
        }                                                                      \
    };                                                                         \
                                                                               \
    template <class T, class U>                                                \
    class __flint_expr<mpq_t,                                                  \
                       __flint_binary_expr<__flint_expr<fmpq_t, T>,            \
                                           __flint_expr<fmpz_t, U>, eval_fun>> \
    {                                                                          \
    private:                                                                   \
        typedef __flint_expr<fmpq_t, T> val1_type;                             \
        typedef __flint_expr<fmpz_t, U> val2_type;                             \
                                                                               \
        __flint_binary_expr<val1_type, val2_type, eval_fun> expr;              \
                                                                               \
    public:                                                                    \
        __flint_expr(const val1_type& val1, const val2_type& val2) :           \
                expr(val1, val2)                                               \
        {}                                                                     \
        void eval(fmpq* q) const                                               \
        {                                                                      \
            fmpz_class temp2(expr.val2);                                       \
            expr.val1.eval(q);                                                 \
            eval_fun::eval(q, q, temp2.get_fmpz_t());                          \
        }                                                                      \
        const val1_type& get_val1() const                                      \
        {                                                                      \
            return expr.val1;                                                  \
        }                                                                      \
        const val2_type& get_val2() const                                      \
        {                                                                      \
            return expr.val2;                                                  \
        }                                                                      \
    };

__FLINTZQ_DEFINE_EXPR(__flint_binary_plus)
__FLINTZQ_DEFINE_EXPR(__flint_binary_minus)

/**************** Macros for defining functions ****************/
/* Results of operators and functions are instances of __flint_expr<T, U>.
   T determines the numerical type of the expression: it can be either
   fmpz_t or fmpq_t.  When the arguments of a binary
   expression have different numerical types, __flint_resolve_expr is used
   to determine the "larger" type.
   U is either __flint_unary_expr<V, Op> or __flint_binary_expr<V, W, Op>,
   where V and W are the arguments' types -- they can in turn be
   expressions, thus allowing to build compound expressions to any
   degree of complexity.
   Op is a function object that must have an eval() method accepting
   appropriate arguments.
   Actual evaluation of a __flint_expr<T, U> object is done when it gets
   assigned to an mp*_class ("lazy" evaluation): this is done by calling
   its eval() method. */

// non-member unary operators and functions

#define __FLINT_DEFINE_UNARY_FUNCTION(fun, eval_fun)                           \
                                                                               \
    template <class T, class U>                                                \
    inline __flint_expr<T, __flint_unary_expr<__flint_expr<T, U>, eval_fun>>   \
    fun(const __flint_expr<T, U>& expr)                                        \
    {                                                                          \
        return __flint_expr<T,                                                 \
                            __flint_unary_expr<__flint_expr<T, U>, eval_fun>>( \
            expr);                                                             \
    }

// variant that only works for one of { mpz, mpq }

#define __FLINT_DEFINE_UNARY_FUNCTION_1(T, fun, eval_fun)                      \
                                                                               \
    template <class U>                                                         \
    inline __flint_expr<T, __flint_unary_expr<__flint_expr<T, U>, eval_fun>>   \
    fun(const __flint_expr<T, U>& expr)                                        \
    {                                                                          \
        return __flint_expr<T,                                                 \
                            __flint_unary_expr<__flint_expr<T, U>, eval_fun>>( \
            expr);                                                             \
    }

#define __FLINT_DEFINE_UNARY_TYPE_FUNCTION(type, fun, eval_fun) \
                                                                \
    template <class T, class U>                                 \
    inline type fun(const __flint_expr<T, U>& expr)             \
    {                                                           \
        __flint_expr<T, T> const& temp(expr);                   \
        return eval_fun::eval(temp.__get_mp());                 \
    }

// non-member binary operators and functions

#define __FLINTP_DEFINE_BINARY_FUNCTION(fun, eval_fun)                         \
                                                                               \
    template <class T, class U, class V, class W>                              \
    inline __flint_expr<                                                       \
        typename __flint_resolve_expr<T, V>::value_type,                       \
        __flint_binary_expr<__flint_expr<T, U>, __flint_expr<V, W>, eval_fun>> \
    fun(const __flint_expr<T, U>& expr1, const __flint_expr<V, W>& expr2)      \
    {                                                                          \
        return __flint_expr<                                                   \
            typename __flint_resolve_expr<T, V>::value_type,                   \
            __flint_binary_expr<__flint_expr<T, U>, __flint_expr<V, W>,        \
                                eval_fun>>(expr1, expr2);                      \
    }

#define __FLINTNN_DEFINE_BINARY_FUNCTION(fun, eval_fun, type, bigtype)      \
                                                                            \
    template <class T, class U>                                             \
    inline __flint_expr<                                                    \
        T, __flint_binary_expr<__flint_expr<T, U>, bigtype, eval_fun>>      \
    fun(const __flint_expr<T, U>& expr, type t)                             \
    {                                                                       \
        return __flint_expr<                                                \
            T, __flint_binary_expr<__flint_expr<T, U>, bigtype, eval_fun>>( \
            expr, t);                                                       \
    }                                                                       \
                                                                            \
    template <class T, class U>                                             \
    inline __flint_expr<                                                    \
        T, __flint_binary_expr<bigtype, __flint_expr<T, U>, eval_fun>>      \
    fun(type t, const __flint_expr<T, U>& expr)                             \
    {                                                                       \
        return __flint_expr<                                                \
            T, __flint_binary_expr<bigtype, __flint_expr<T, U>, eval_fun>>( \
            t, expr);                                                       \
    }

#define __FLINTNS_DEFINE_BINARY_FUNCTION(fun, eval_fun, type) \
    __FLINTNN_DEFINE_BINARY_FUNCTION(fun, eval_fun, type, signed long int)

#define __FLINTNU_DEFINE_BINARY_FUNCTION(fun, eval_fun, type) \
    __FLINTNN_DEFINE_BINARY_FUNCTION(fun, eval_fun, type, unsigned long int)

#define __FLINTND_DEFINE_BINARY_FUNCTION(fun, eval_fun, type) \
    __FLINTNN_DEFINE_BINARY_FUNCTION(fun, eval_fun, type, double)

#define __FLINTNLD_DEFINE_BINARY_FUNCTION(fun, eval_fun, type) \
    __FLINTNN_DEFINE_BINARY_FUNCTION(fun, eval_fun, type, long double)

#define __FLINTN_DEFINE_BINARY_FUNCTION(fun, eval_fun)                  \
    __FLINTNS_DEFINE_BINARY_FUNCTION(fun, eval_fun, signed char)        \
    __FLINTNU_DEFINE_BINARY_FUNCTION(fun, eval_fun, unsigned char)      \
    __FLINTNS_DEFINE_BINARY_FUNCTION(fun, eval_fun, signed int)         \
    __FLINTNU_DEFINE_BINARY_FUNCTION(fun, eval_fun, unsigned int)       \
    __FLINTNS_DEFINE_BINARY_FUNCTION(fun, eval_fun, signed short int)   \
    __FLINTNU_DEFINE_BINARY_FUNCTION(fun, eval_fun, unsigned short int) \
    __FLINTNS_DEFINE_BINARY_FUNCTION(fun, eval_fun, signed long int)    \
    __FLINTNU_DEFINE_BINARY_FUNCTION(fun, eval_fun, unsigned long int)  \
    __FLINTND_DEFINE_BINARY_FUNCTION(fun, eval_fun, float)              \
    __FLINTND_DEFINE_BINARY_FUNCTION(fun, eval_fun, double)             \
    /* __FLINTNLD_DEFINE_BINARY_FUNCTION(fun, eval_fun, long double) */

#define __FLINT_DEFINE_BINARY_FUNCTION(fun, eval_fun) \
    __FLINTP_DEFINE_BINARY_FUNCTION(fun, eval_fun)    \
    __FLINTN_DEFINE_BINARY_FUNCTION(fun, eval_fun)

// variant that only works for one of { mpz, mpq }

#define __FLINTP_DEFINE_BINARY_FUNCTION_1(T, fun, eval_fun)                   \
                                                                              \
    template <class U, class W>                                               \
    inline __flint_expr<T, __flint_binary_expr<__flint_expr<T, U>,            \
                                               __flint_expr<T, W>, eval_fun>> \
    fun(const __flint_expr<T, U>& expr1, const __flint_expr<T, W>& expr2)     \
    {                                                                         \
        return __flint_expr<                                                  \
            T, __flint_binary_expr<__flint_expr<T, U>, __flint_expr<T, W>,    \
                                   eval_fun>>(expr1, expr2);                  \
    }

#define __FLINTNN_DEFINE_BINARY_FUNCTION_1(T, fun, eval_fun, type, bigtype) \
                                                                            \
    template <class U>                                                      \
    inline __flint_expr<                                                    \
        T, __flint_binary_expr<__flint_expr<T, U>, bigtype, eval_fun>>      \
    fun(const __flint_expr<T, U>& expr, type t)                             \
    {                                                                       \
        return __flint_expr<                                                \
            T, __flint_binary_expr<__flint_expr<T, U>, bigtype, eval_fun>>( \
            expr, t);                                                       \
    }                                                                       \
                                                                            \
    template <class U>                                                      \
    inline __flint_expr<                                                    \
        T, __flint_binary_expr<bigtype, __flint_expr<T, U>, eval_fun>>      \
    fun(type t, const __flint_expr<T, U>& expr)                             \
    {                                                                       \
        return __flint_expr<                                                \
            T, __flint_binary_expr<bigtype, __flint_expr<T, U>, eval_fun>>( \
            t, expr);                                                       \
    }

#define __FLINTNS_DEFINE_BINARY_FUNCTION_1(T, fun, eval_fun, type) \
    __FLINTNN_DEFINE_BINARY_FUNCTION_1(T, fun, eval_fun, type, signed long int)

#define __FLINTNU_DEFINE_BINARY_FUNCTION_1(T, fun, eval_fun, type) \
    __FLINTNN_DEFINE_BINARY_FUNCTION_1(T, fun, eval_fun, type,     \
                                       unsigned long int)

#define __FLINTND_DEFINE_BINARY_FUNCTION_1(T, fun, eval_fun, type) \
    __FLINTNN_DEFINE_BINARY_FUNCTION_1(T, fun, eval_fun, type, double)

#define __FLINTNLD_DEFINE_BINARY_FUNCTION_1(T, fun, eval_fun, type) \
    __FLINTNN_DEFINE_BINARY_FUNCTION_1(T, fun, eval_fun, type, long double)

#define __FLINTN_DEFINE_BINARY_FUNCTION_1(T, fun, eval_fun)                  \
    __FLINTNS_DEFINE_BINARY_FUNCTION_1(T, fun, eval_fun, signed char)        \
    __FLINTNU_DEFINE_BINARY_FUNCTION_1(T, fun, eval_fun, unsigned char)      \
    __FLINTNS_DEFINE_BINARY_FUNCTION_1(T, fun, eval_fun, signed int)         \
    __FLINTNU_DEFINE_BINARY_FUNCTION_1(T, fun, eval_fun, unsigned int)       \
    __FLINTNS_DEFINE_BINARY_FUNCTION_1(T, fun, eval_fun, signed short int)   \
    __FLINTNU_DEFINE_BINARY_FUNCTION_1(T, fun, eval_fun, unsigned short int) \
    __FLINTNS_DEFINE_BINARY_FUNCTION_1(T, fun, eval_fun, signed long int)    \
    __FLINTNU_DEFINE_BINARY_FUNCTION_1(T, fun, eval_fun, unsigned long int)  \
    __FLINTND_DEFINE_BINARY_FUNCTION_1(T, fun, eval_fun, float)              \
    __FLINTND_DEFINE_BINARY_FUNCTION_1(T, fun, eval_fun, double)

#define __FLINT_DEFINE_BINARY_FUNCTION_1(T, fun, eval_fun) \
    __FLINTP_DEFINE_BINARY_FUNCTION_1(T, fun, eval_fun)    \
    __FLINTN_DEFINE_BINARY_FUNCTION_1(T, fun, eval_fun)

#define __FLINT_DEFINE_BINARY_FUNCTION_UI(fun, eval_fun)                    \
                                                                            \
    template <class T, class U>                                             \
    inline __flint_expr<                                                    \
        T, __flint_binary_expr<__flint_expr<T, U>, mp_bitcnt_t, eval_fun>>  \
    fun(const __flint_expr<T, U>& expr, mp_bitcnt_t l)                      \
    {                                                                       \
        return __flint_expr<T, __flint_binary_expr<__flint_expr<T, U>,      \
                                                   mp_bitcnt_t, eval_fun>>( \
            expr, l);                                                       \
    }

#define __FLINTP_DEFINE_BINARY_TYPE_FUNCTION(type, fun, eval_fun)  \
                                                                   \
    template <class T, class U, class V, class W>                  \
    inline type fun(const __flint_expr<T, U>& expr1,               \
                    const __flint_expr<V, W>& expr2)               \
    {                                                              \
        __flint_expr<T, T> const& temp1(expr1);                    \
        __flint_expr<V, V> const& temp2(expr2);                    \
        return eval_fun::eval(temp1.__get_mp(), temp2.__get_mp()); \
    }

#define __FLINTNN_DEFINE_BINARY_TYPE_FUNCTION(type, fun, eval_fun, type2, \
                                              bigtype)                    \
                                                                          \
    template <class T, class U>                                           \
    inline type fun(const __flint_expr<T, U>& expr, type2 t)              \
    {                                                                     \
        __flint_expr<T, T> const& temp(expr);                             \
        return eval_fun::eval(temp.__get_mp(), static_cast<bigtype>(t));  \
    }                                                                     \
                                                                          \
    template <class T, class U>                                           \
    inline type fun(type2 t, const __flint_expr<T, U>& expr)              \
    {                                                                     \
        __flint_expr<T, T> const& temp(expr);                             \
        return eval_fun::eval(static_cast<bigtype>(t), temp.__get_mp());  \
    }

#define __FLINTNS_DEFINE_BINARY_TYPE_FUNCTION(type, fun, eval_fun, type2) \
    __FLINTNN_DEFINE_BINARY_TYPE_FUNCTION(type, fun, eval_fun, type2,     \
                                          signed long int)

#define __FLINTNU_DEFINE_BINARY_TYPE_FUNCTION(type, fun, eval_fun, type2) \
    __FLINTNN_DEFINE_BINARY_TYPE_FUNCTION(type, fun, eval_fun, type2,     \
                                          unsigned long int)

#define __FLINTND_DEFINE_BINARY_TYPE_FUNCTION(type, fun, eval_fun, type2) \
    __FLINTNN_DEFINE_BINARY_TYPE_FUNCTION(type, fun, eval_fun, type2, double)

#define __FLINTNLD_DEFINE_BINARY_TYPE_FUNCTION(type, fun, eval_fun, type2) \
    __FLINTNN_DEFINE_BINARY_TYPE_FUNCTION(type, fun, eval_fun, type2,      \
                                          long double)

#define __FLINTN_DEFINE_BINARY_TYPE_FUNCTION(type, fun, eval_fun)             \
    __FLINTNS_DEFINE_BINARY_TYPE_FUNCTION(type, fun, eval_fun, signed char)   \
    __FLINTNU_DEFINE_BINARY_TYPE_FUNCTION(type, fun, eval_fun, unsigned char) \
    __FLINTNS_DEFINE_BINARY_TYPE_FUNCTION(type, fun, eval_fun, signed int)    \
    __FLINTNU_DEFINE_BINARY_TYPE_FUNCTION(type, fun, eval_fun, unsigned int)  \
    __FLINTNS_DEFINE_BINARY_TYPE_FUNCTION(type, fun, eval_fun,                \
                                          signed short int)                   \
    __FLINTNU_DEFINE_BINARY_TYPE_FUNCTION(type, fun, eval_fun,                \
                                          unsigned short int)                 \
    __FLINTNS_DEFINE_BINARY_TYPE_FUNCTION(type, fun, eval_fun,                \
                                          signed long int)                    \
    __FLINTNU_DEFINE_BINARY_TYPE_FUNCTION(type, fun, eval_fun,                \
                                          unsigned long int)                  \
    __FLINTND_DEFINE_BINARY_TYPE_FUNCTION(type, fun, eval_fun, float)         \
    __FLINTND_DEFINE_BINARY_TYPE_FUNCTION(type, fun, eval_fun, double)

#define __FLINT_DEFINE_BINARY_TYPE_FUNCTION(type, fun, eval_fun) \
    __FLINTP_DEFINE_BINARY_TYPE_FUNCTION(type, fun, eval_fun)    \
    __FLINTN_DEFINE_BINARY_TYPE_FUNCTION(type, fun, eval_fun)

// member operators

#define __FLINTP_DEFINE_COMPOUND_OPERATOR(type, fun, eval_fun)                 \
                                                                               \
    template <class T, class U>                                                \
    inline type##_class& type##_class::fun(__flint_expr<T, U> const& expr)     \
    {                                                                          \
        __flint_set_expr(                                                      \
            &mp,                                                               \
            __flint_expr<type##_t,                                             \
                         __flint_binary_expr<type##_class, __flint_expr<T, U>, \
                                             eval_fun>>(*this, expr));         \
        return *this;                                                          \
    }

#define __FLINTNN_DEFINE_COMPOUND_OPERATOR(type, fun, eval_fun, type2,        \
                                           bigtype)                           \
                                                                              \
    inline type##_class& type##_class::fun(type2 t)                           \
    {                                                                         \
        __flint_set_expr(                                                     \
            &mp,                                                              \
            __flint_expr<type##_t, __flint_binary_expr<type##_class, bigtype, \
                                                       eval_fun>>(*this, t)); \
        return *this;                                                         \
    }

#define __FLINTNS_DEFINE_COMPOUND_OPERATOR(type, fun, eval_fun, type2) \
    __FLINTNN_DEFINE_COMPOUND_OPERATOR(type, fun, eval_fun, type2,     \
                                       signed long int)

#define __FLINTNU_DEFINE_COMPOUND_OPERATOR(type, fun, eval_fun, type2) \
    __FLINTNN_DEFINE_COMPOUND_OPERATOR(type, fun, eval_fun, type2,     \
                                       unsigned long int)

#define __FLINTND_DEFINE_COMPOUND_OPERATOR(type, fun, eval_fun, type2) \
    __FLINTNN_DEFINE_COMPOUND_OPERATOR(type, fun, eval_fun, type2, double)

#define __FLINTNLD_DEFINE_COMPOUND_OPERATOR(type, fun, eval_fun, type2) \
    __FLINTNN_DEFINE_COMPOUND_OPERATOR(type, fun, eval_fun, type2, long double)

#define __FLINTN_DEFINE_COMPOUND_OPERATOR(type, fun, eval_fun)                 \
    __FLINTNS_DEFINE_COMPOUND_OPERATOR(type, fun, eval_fun, signed char)       \
    __FLINTNU_DEFINE_COMPOUND_OPERATOR(type, fun, eval_fun, unsigned char)     \
    __FLINTNS_DEFINE_COMPOUND_OPERATOR(type, fun, eval_fun, signed int)        \
    __FLINTNU_DEFINE_COMPOUND_OPERATOR(type, fun, eval_fun, unsigned int)      \
    __FLINTNS_DEFINE_COMPOUND_OPERATOR(type, fun, eval_fun, signed short int)  \
    __FLINTNU_DEFINE_COMPOUND_OPERATOR(type, fun, eval_fun,                    \
                                       unsigned short int)                     \
    __FLINTNS_DEFINE_COMPOUND_OPERATOR(type, fun, eval_fun, signed long int)   \
    __FLINTNU_DEFINE_COMPOUND_OPERATOR(type, fun, eval_fun, unsigned long int) \
    __FLINTND_DEFINE_COMPOUND_OPERATOR(type, fun, eval_fun, float)             \
    __FLINTND_DEFINE_COMPOUND_OPERATOR(type, fun, eval_fun, double)

#define __FLINT_DEFINE_COMPOUND_OPERATOR(type, fun, eval_fun) \
    __FLINTP_DEFINE_COMPOUND_OPERATOR(type, fun, eval_fun)    \
    __FLINTN_DEFINE_COMPOUND_OPERATOR(type, fun, eval_fun)

#define __FLINTZ_DEFINE_COMPOUND_OPERATOR(fun, eval_fun) \
    __FLINT_DEFINE_COMPOUND_OPERATOR(fmpz, fun, eval_fun)

#define __FLINTQ_DEFINE_COMPOUND_OPERATOR(fun, eval_fun) \
    __FLINT_DEFINE_COMPOUND_OPERATOR(fmpq, fun, eval_fun)

#define __FLINT_DEFINE_COMPOUND_OPERATOR_UI(type, fun, eval_fun)              \
                                                                              \
    inline type##_class& type##_class::fun(mp_bitcnt_t l)                     \
    {                                                                         \
        __flint_set_expr(                                                     \
            &mp,                                                              \
            __flint_expr<type##_t, __flint_binary_expr<                       \
                                       type##_class, mp_bitcnt_t, eval_fun>>( \
                *this, l));                                                   \
        return *this;                                                         \
    }

#define __FLINTZ_DEFINE_COMPOUND_OPERATOR_UI(fun, eval_fun) \
    __FLINT_DEFINE_COMPOUND_OPERATOR_UI(fmpz, fun, eval_fun)

#define __FLINTQ_DEFINE_COMPOUND_OPERATOR_UI(fun, eval_fun) \
    __FLINT_DEFINE_COMPOUND_OPERATOR_UI(fmpq, fun, eval_fun)

#define __FLINT_DEFINE_INCREMENT_OPERATOR(type, fun, eval_fun) \
                                                               \
    inline type##_class& type##_class::fun()                   \
    {                                                          \
        eval_fun::eval(&mp);                                   \
        return *this;                                          \
    }                                                          \
                                                               \
    inline type##_class type##_class::fun(int)                 \
    {                                                          \
        type##_class temp(*this);                              \
        eval_fun::eval(&mp);                                   \
        return temp;                                           \
    }

#define __FLINTZ_DEFINE_INCREMENT_OPERATOR(fun, eval_fun) \
    __FLINT_DEFINE_INCREMENT_OPERATOR(fmpz, fun, eval_fun)

#define __FLINTQ_DEFINE_INCREMENT_OPERATOR(fun, eval_fun) \
    __FLINT_DEFINE_INCREMENT_OPERATOR(fmpq, fun, eval_fun)

#define __FLINTP_DEFINE_UNARY_STATIC_MEMFUN(T, fun, eval_fun)                  \
    template <class U>                                                         \
    __flint_expr<T, __flint_unary_expr<__flint_expr<T, U>, eval_fun>> fun(     \
        const __flint_expr<T, U>& expr)                                        \
    {                                                                          \
        return __flint_expr<T,                                                 \
                            __flint_unary_expr<__flint_expr<T, U>, eval_fun>>( \
            expr);                                                             \
    }

#define __FLINTNN_DEFINE_UNARY_STATIC_MEMFUN(T, fun, eval_fun, type, bigtype) \
    inline __flint_expr<T, __flint_unary_expr<bigtype, eval_fun>> fun(        \
        type expr)                                                            \
    {                                                                         \
        return __flint_expr<T, __flint_unary_expr<bigtype, eval_fun>>(expr);  \
    }

#define __FLINTNS_DEFINE_UNARY_STATIC_MEMFUN(T, fun, eval_fun, type) \
    __FLINTNN_DEFINE_UNARY_STATIC_MEMFUN(T, fun, eval_fun, type, signed long)
#define __FLINTNU_DEFINE_UNARY_STATIC_MEMFUN(T, fun, eval_fun, type) \
    __FLINTNN_DEFINE_UNARY_STATIC_MEMFUN(T, fun, eval_fun, type, unsigned long)
#define __FLINTND_DEFINE_UNARY_STATIC_MEMFUN(T, fun, eval_fun, type) \
    __FLINTNN_DEFINE_UNARY_STATIC_MEMFUN(T, fun, eval_fun, type, double)

#define __FLINTN_DEFINE_UNARY_STATIC_MEMFUN(T, fun, eval_fun)                  \
    __FLINTNS_DEFINE_UNARY_STATIC_MEMFUN(T, fun, eval_fun, signed char)        \
    __FLINTNU_DEFINE_UNARY_STATIC_MEMFUN(T, fun, eval_fun, unsigned char)      \
    __FLINTNS_DEFINE_UNARY_STATIC_MEMFUN(T, fun, eval_fun, signed int)         \
    __FLINTNU_DEFINE_UNARY_STATIC_MEMFUN(T, fun, eval_fun, unsigned int)       \
    __FLINTNS_DEFINE_UNARY_STATIC_MEMFUN(T, fun, eval_fun, signed short int)   \
    __FLINTNU_DEFINE_UNARY_STATIC_MEMFUN(T, fun, eval_fun, unsigned short int) \
    __FLINTNS_DEFINE_UNARY_STATIC_MEMFUN(T, fun, eval_fun, signed long int)    \
    __FLINTNU_DEFINE_UNARY_STATIC_MEMFUN(T, fun, eval_fun, unsigned long int)  \
    __FLINTND_DEFINE_UNARY_STATIC_MEMFUN(T, fun, eval_fun, float)              \
    __FLINTND_DEFINE_UNARY_STATIC_MEMFUN(T, fun, eval_fun, double)

#define __FLINT_DEFINE_UNARY_STATIC_MEMFUN(T, fun, eval_fun) \
    __FLINTP_DEFINE_UNARY_STATIC_MEMFUN(T, fun, eval_fun)    \
    __FLINTN_DEFINE_UNARY_STATIC_MEMFUN(T, fun, eval_fun)

/**************** Arithmetic operators and functions ****************/

// non-member operators and functions

__FLINT_DEFINE_UNARY_FUNCTION(operator+, __flint_unary_plus)

__FLINT_DEFINE_UNARY_FUNCTION(operator-, __flint_unary_minus)
__FLINT_DEFINE_UNARY_FUNCTION_1(fmpz_t, operator~, __flint_unary_com)

__FLINT_DEFINE_BINARY_FUNCTION(operator+, __flint_binary_plus)

__FLINT_DEFINE_BINARY_FUNCTION(operator-, __flint_binary_minus)

__FLINT_DEFINE_BINARY_FUNCTION(operator*, __flint_binary_multiplies)

__FLINT_DEFINE_BINARY_FUNCTION(operator/, __flint_binary_divides)

__FLINT_DEFINE_BINARY_FUNCTION_1(fmpz_t, operator%, __flint_binary_modulus)

__FLINT_DEFINE_BINARY_FUNCTION_1(fmpz_t, operator&, __flint_binary_and)

__FLINT_DEFINE_BINARY_FUNCTION_1(fmpz_t, operator|, __flint_binary_ior)
__FLINT_DEFINE_BINARY_FUNCTION_1(fmpz_t, operator^, __flint_binary_xor)

__FLINT_DEFINE_BINARY_FUNCTION_UI(operator<<, __flint_binary_lshift)
__FLINT_DEFINE_BINARY_FUNCTION_UI(operator>>, __flint_binary_rshift)

__FLINT_DEFINE_BINARY_TYPE_FUNCTION(bool, operator==, __flint_binary_equal)

__FLINT_DEFINE_BINARY_TYPE_FUNCTION(bool, operator!=, !__flint_binary_equal)

__FLINT_DEFINE_BINARY_TYPE_FUNCTION(bool, operator<, __flint_binary_less)

__FLINT_DEFINE_BINARY_TYPE_FUNCTION(bool, operator<=, !__flint_binary_greater)

__FLINT_DEFINE_BINARY_TYPE_FUNCTION(bool, operator>, __flint_binary_greater)
__FLINT_DEFINE_BINARY_TYPE_FUNCTION(bool, operator>=, !__flint_binary_less)

__FLINT_DEFINE_UNARY_FUNCTION(abs, __flint_abs_function)

__FLINT_DEFINE_UNARY_FUNCTION_1(fmpz_t, sqrt, __flint_sqrt_function)
__FLINT_DEFINE_BINARY_FUNCTION_1(fmpz_t, gcd, __flint_gcd_function)
__FLINT_DEFINE_BINARY_FUNCTION_1(fmpz_t, lcm, __flint_lcm_function)

__FLINT_DEFINE_UNARY_TYPE_FUNCTION(int, sgn, __flint_sgn_function)
__FLINT_DEFINE_BINARY_TYPE_FUNCTION(int, cmp, __flint_cmp_function)

template <class T>
void swap(__flint_expr<T, T>& x, __flint_expr<T, T>& y) noexcept
{
    x.swap(y);
}

// member operators for fmpz_class

__FLINTZ_DEFINE_COMPOUND_OPERATOR(operator+=, __flint_binary_plus)
__FLINTZ_DEFINE_COMPOUND_OPERATOR(operator-=, __flint_binary_minus)
__FLINTZ_DEFINE_COMPOUND_OPERATOR(operator*=, __flint_binary_multiplies)
__FLINTZ_DEFINE_COMPOUND_OPERATOR(operator/=, __flint_binary_divides)
__FLINTZ_DEFINE_COMPOUND_OPERATOR(operator%=, __flint_binary_modulus)

__FLINTZ_DEFINE_COMPOUND_OPERATOR(operator&=, __flint_binary_and)
__FLINTZ_DEFINE_COMPOUND_OPERATOR(operator|=, __flint_binary_ior)
__FLINTZ_DEFINE_COMPOUND_OPERATOR(operator^=, __flint_binary_xor)

__FLINTZ_DEFINE_COMPOUND_OPERATOR_UI(operator<<=, __flint_binary_lshift)
__FLINTZ_DEFINE_COMPOUND_OPERATOR_UI(operator>>=, __flint_binary_rshift)

__FLINTZ_DEFINE_INCREMENT_OPERATOR(operator++, __flint_unary_increment)
__FLINTZ_DEFINE_INCREMENT_OPERATOR(operator--, __flint_unary_decrement)

// member operators for fmpq_class

__FLINTQ_DEFINE_COMPOUND_OPERATOR(operator+=, __flint_binary_plus)
__FLINTQ_DEFINE_COMPOUND_OPERATOR(operator-=, __flint_binary_minus)
__FLINTQ_DEFINE_COMPOUND_OPERATOR(operator*=, __flint_binary_multiplies)
__FLINTQ_DEFINE_COMPOUND_OPERATOR(operator/=, __flint_binary_divides)

__FLINTQ_DEFINE_COMPOUND_OPERATOR_UI(operator<<=, __flint_binary_lshift)
__FLINTQ_DEFINE_COMPOUND_OPERATOR_UI(operator>>=, __flint_binary_rshift)

__FLINTQ_DEFINE_INCREMENT_OPERATOR(operator++, __flint_unary_increment)
__FLINTQ_DEFINE_INCREMENT_OPERATOR(operator--, __flint_unary_decrement)

/**************** Specialize std::numeric_limits ****************/

namespace std
{

template <>
class numeric_limits<fmpz_class>
{
public:
    static bool const is_specialized = true;

    static fmpz_class min() { return {}; }

    static fmpz_class max() { return {}; }

    static fmpz_class lowest() { return {}; }

    static int const digits       = 0;
    static int const digits10     = 0;
    static int const max_digits10 = 0;
    static bool const is_signed   = true;
    static bool const is_integer  = true;
    static bool const is_exact    = true;
    static int const radix        = 2;

    static fmpz_class epsilon() { return {}; }

    static fmpz_class round_error() { return {}; }

    static int const min_exponent              = 0;
    static int const min_exponent10            = 0;
    static int const max_exponent              = 0;
    static int const max_exponent10            = 0;
    static bool const has_infinity             = false;
    static bool const has_quiet_NaN            = false;
    static bool const has_signaling_NaN        = false;
    static float_denorm_style const has_denorm = denorm_absent;
    static bool const has_denorm_loss          = false;

    static fmpz_class infinity() { return {}; }

    static fmpz_class quiet_NaN() { return {}; }

    static fmpz_class signaling_NaN() { return {}; }

    static fmpz_class denorm_min() { return {}; }

    static bool const is_iec559                = false;
    static bool const is_bounded               = false;
    static bool const is_modulo                = false;
    static bool const traps                    = false;
    static bool const tinyness_before          = false;
    static float_round_style const round_style = round_toward_zero;
};

template <>
class numeric_limits<fmpq_class>
{
public:
    static bool const is_specialized = true;

    static fmpq_class min() { return {}; }

    static fmpq_class max() { return {}; }

    static fmpq_class lowest() { return {}; }

    static int const digits       = 0;
    static int const digits10     = 0;
    static int const max_digits10 = 0;
    static bool const is_signed   = true;
    static bool const is_integer  = false;
    static bool const is_exact    = true;
    static int const radix        = 2;

    static fmpq_class epsilon() { return {}; }

    static fmpq_class round_error() { return {}; }

    static int const min_exponent              = 0;
    static int const min_exponent10            = 0;
    static int const max_exponent              = 0;
    static int const max_exponent10            = 0;
    static bool const has_infinity             = false;
    static bool const has_quiet_NaN            = false;
    static bool const has_signaling_NaN        = false;
    static float_denorm_style const has_denorm = denorm_absent;
    static bool const has_denorm_loss          = false;

    static fmpq_class infinity() { return {}; }

    static fmpq_class quiet_NaN() { return {}; }

    static fmpq_class signaling_NaN() { return {}; }

    static fmpq_class denorm_min() { return {}; }

    static bool const is_iec559                = false;
    static bool const is_bounded               = false;
    static bool const is_modulo                = false;
    static bool const traps                    = false;
    static bool const tinyness_before          = false;
    static float_round_style const round_style = round_toward_zero;
};

}  // namespace std

/**************** #undef all private macros ****************/

#undef __FLINTP_DECLARE_COMPOUND_OPERATOR
#undef __FLINTN_DECLARE_COMPOUND_OPERATOR
#undef __FLINT_DECLARE_COMPOUND_OPERATOR
#undef __FLINT_DECLARE_COMPOUND_OPERATOR_UI
#undef __FLINT_DECLARE_INCREMENT_OPERATOR
#undef __FLINTXX_DEFINE_ARITHMETIC_CONSTRUCTORS
#undef __FLINTXX_DEFINE_ARITHMETIC_ASSIGNMENTS

#undef __FLINTZQ_DEFINE_EXPR

#undef __FLINT_DEFINE_UNARY_FUNCTION_1
#undef __FLINT_DEFINE_UNARY_FUNCTION
#undef __FLINT_DEFINE_UNARY_TYPE_FUNCTION

#undef __FLINTP_DEFINE_BINARY_FUNCTION
#undef __FLINTNN_DEFINE_BINARY_FUNCTION
#undef __FLINTNS_DEFINE_BINARY_FUNCTION
#undef __FLINTNU_DEFINE_BINARY_FUNCTION
#undef __FLINTND_DEFINE_BINARY_FUNCTION
#undef __FLINTNLD_DEFINE_BINARY_FUNCTION
#undef __FLINTN_DEFINE_BINARY_FUNCTION
#undef __FLINT_DEFINE_BINARY_FUNCTION

#undef __FLINT_DEFINE_BINARY_FUNCTION_UI

#undef __FLINTP_DEFINE_BINARY_TYPE_FUNCTION
#undef __FLINTNN_DEFINE_BINARY_TYPE_FUNCTION
#undef __FLINTNS_DEFINE_BINARY_TYPE_FUNCTION
#undef __FLINTNU_DEFINE_BINARY_TYPE_FUNCTION
#undef __FLINTND_DEFINE_BINARY_TYPE_FUNCTION
#undef __FLINTNLD_DEFINE_BINARY_TYPE_FUNCTION
#undef __FLINTN_DEFINE_BINARY_TYPE_FUNCTION
#undef __FLINT_DEFINE_BINARY_TYPE_FUNCTION

#undef __FLINTZ_DEFINE_COMPOUND_OPERATOR

#undef __FLINTP_DEFINE_COMPOUND_OPERATOR
#undef __FLINTNN_DEFINE_COMPOUND_OPERATOR
#undef __FLINTNS_DEFINE_COMPOUND_OPERATOR
#undef __FLINTNU_DEFINE_COMPOUND_OPERATOR
#undef __FLINTND_DEFINE_COMPOUND_OPERATOR
#undef __FLINTNLD_DEFINE_COMPOUND_OPERATOR
#undef __FLINTN_DEFINE_COMPOUND_OPERATOR
#undef __FLINT_DEFINE_COMPOUND_OPERATOR

#undef __FLINTQ_DEFINE_COMPOUND_OPERATOR
#undef __FLINTF_DEFINE_COMPOUND_OPERATOR

#undef __FLINT_DEFINE_COMPOUND_OPERATOR_UI
#undef __FLINTZ_DEFINE_COMPOUND_OPERATOR_UI
#undef __FLINTQ_DEFINE_COMPOUND_OPERATOR_UI
#undef __FLINTF_DEFINE_COMPOUND_OPERATOR_UI

#undef __FLINT_DEFINE_INCREMENT_OPERATOR
#undef __FLINTZ_DEFINE_INCREMENT_OPERATOR
#undef __FLINTQ_DEFINE_INCREMENT_OPERATOR
#undef __FLINTF_DEFINE_INCREMENT_OPERATOR

#undef __FLINTXX_CONSTANT_TRUE
#undef __FLINTXX_CONSTANT

// NOLINTEND(cert-dcl51-cpp,cert-dcl37-c,bugprone-reserved-identifier)
// NOLINTEND(google-explicit-constructor,hicpp-explicit-conversions)
// NOLINTEND(bugprone-macro-parentheses)
// NOLINTEND(google-runtime-int)
