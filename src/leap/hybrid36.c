// Derived from:
// C port of hy36encode/hy36decode. See hybrid_36.py for details.
// cctbx.sf.net — Ralf W. Grosse-Kunstleve, Feb 2007.
#include "hybrid36.h"


// NOTE: full printf format would support sign and zero fill.
// Here we use zero fill for overflow wrapped integers.
/* ------------------------------------------------------------------ */
/* Module-level constants                                               */
/* ------------------------------------------------------------------ */

static const char DIGITS_UPPER[] = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";
static const char DIGITS_LOWER[] = "0123456789abcdefghijklmnopqrstuvwxyz";

/* ------------------------------------------------------------------ */
/* Internal helpers                                                     */
/* ------------------------------------------------------------------ */

static void fill_with_asterisks(unsigned width, char *result)
{
    while (width--) *result++ = '*';
    *result = '\0';
}

static void encode_pure(
    const char *digits,
    unsigned    digits_size,
    unsigned    width,
    int         value,
    char        fill_char,
    char       *result)
{
    char     buf[16];
    unsigned i = 0;
    unsigned pad;
    int      negative = (value < 0);

    if (negative) value = -value;

    do {
        int rest = value / (int)digits_size;
        buf[i++] = digits[value - rest * (int)digits_size];
        value = rest;
    } while (value != 0);

    if (negative) buf[i++] = '-';

    for (pad = i; pad < width; pad++) *result++ = fill_char;
    while (i != 0) *result++ = buf[--i];
    *result = '\0';
}

static Hy36Error decode_pure(
    const int  *digits_values,
    unsigned    digits_size,
    const char *s,
    unsigned    s_size,
    int        *result)
{
    int      value        = 0;
    int      have_minus   = 0;
    int      have_nonblank = 0;
    unsigned i;

    for (i = 0; i < s_size; i++) {
        int si = (unsigned char)s[i];
        if (si > 127) { *result = 0; return HY36_INVALID_LITERAL; }

        // XXX This allows nonleading blanks to count as zero. Probably a bad idea.
        if (si == ' ') {
            if (have_nonblank) value *= (int)digits_size;
        } else if (si == '-') {
            if (have_nonblank) { *result = 0; return HY36_INVALID_LITERAL; }
            have_nonblank = 1;
            have_minus    = 1;
        } else {
            int dv = digits_values[si];
            if (dv < 0 || dv >= (int)digits_size) {
                *result = 0;
                return HY36_INVALID_LITERAL;
            }
            have_nonblank = 1;
            value = value * (int)digits_size + dv;
        }
    }

    *result = have_minus ? -value : value;
    return HY36_OK;
}

/* ------------------------------------------------------------------ */
/* Lookup table initialisation (called once)                            */
/* ------------------------------------------------------------------ */

static int digits_values_upper[128];
static int digits_values_lower[128];

static Hy36Error init_decode_tables(void)
{
    unsigned i;
    for (i = 0; i < 128U; i++) digits_values_upper[i] = digits_values_lower[i] = -1;

    for (i = 0; i < 36U; i++) {
        int du = (unsigned char)DIGITS_UPPER[i];
        int dl = (unsigned char)DIGITS_LOWER[i];
        if (du > 127 || dl > 127) return HY36_INTERNAL_ERROR;
        digits_values_upper[du] = (int)i;
        digits_values_lower[dl] = (int)i;
    }
    return HY36_OK;
}

/* ------------------------------------------------------------------ */
/* Public API                                                           */
/* ------------------------------------------------------------------ */
Hy36Error hy36encode(unsigned width, int value, char *result)
{
    int i = value;
    int status = HY36_OK;
    char fill_char = ' ';

    if (width == 4U) {
        if (i >= -999) {
            while (i>HY36_WIDTH_4_MAX) {
                i -= (HY36_WIDTH_4_MAX+1);
                // wrap overflow back to zero, use '0' fill char to hint at trucnation
                // Output is valid format but truncated in value
                fill_char = '0';
                status = HY36_WRAPPED;
            }
            if (i < 10000) {
                encode_pure(DIGITS_UPPER, 10U, 4U, i, fill_char, result);
                return status;
            }
            i -= 10000;
            if (i < 1213056 /* 26*36^3 */) {
                encode_pure(DIGITS_UPPER, 36U, 0U, i + 466560 /* 10*36^3 */, fill_char, result);
                return status;
            }
            i -= 1213056;
            encode_pure(DIGITS_LOWER, 36U, 0U, i + 466560, fill_char, result);
            return status;
        }
    }
 else if (width == 5U) {
        if (i >= -9999) {
            while (i>HY36_WIDTH_5_MAX) {
                i -= (HY36_WIDTH_5_MAX+1);
                fill_char = '0';
                status = HY36_WRAPPED;
            }
            if (i < 100000) {
                encode_pure(DIGITS_UPPER, 10U, 5U, i, fill_char, result);
                return status;
            }
            i -= 100000;
            if (i < 43670016 /* 26*36^4 */) {
                encode_pure(DIGITS_UPPER, 36U, 0U, i + 16796160 /* 10*36^4 */, fill_char, result);
                return status;
            }
            i -= 43670016;
            encode_pure(DIGITS_LOWER, 36U, 0U, i + 16796160, fill_char, result);
            return status;
        }
    } else {
        fill_with_asterisks(width, result);
        return HY36_UNSUPPORTED_WIDTH;
    }
    fill_with_asterisks(width, result);
    return HY36_OUT_OF_RANGE;
}

Hy36Error hy36decode(unsigned width, const char *s, int *result)
{
    static int tables_ready = 0;
    int        di;
    Hy36Error  err;

    if (!tables_ready) {
        err = init_decode_tables();
        if (err) { *result = 0; return err; }
        tables_ready = 1;
    }

    di = (unsigned char)s[0];
    if (di > 127)          { *result = 0; return HY36_INVALID_LITERAL; }

    if (digits_values_upper[di] >= 10) {
        err = decode_pure(digits_values_upper, 36U, s, width, result);
        if (err) return err;
        if      (width == 4U) *result -= 456560;
        else if (width == 5U) *result -= 16696160;
        else { *result = 0; return HY36_UNSUPPORTED_WIDTH; }
        return HY36_OK;
    }

    if (digits_values_lower[di] >= 10) {
        err = decode_pure(digits_values_lower, 36U, s, width, result);
        if (err) return err;
        if      (width == 4U) *result += 756496;
        else if (width == 5U) *result += 26973856;
        else { *result = 0; return HY36_UNSUPPORTED_WIDTH; }
        return HY36_OK;
    }

    err = decode_pure(digits_values_upper, 10U, s, width, result);
    if (err) return err;
    if (width != 4U && width != 5U) { *result = 0; return HY36_UNSUPPORTED_WIDTH; }
    return HY36_OK;
}

