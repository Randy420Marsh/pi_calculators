#include <gmp.h>
#include <mpfr.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <limits.h>

/* =========================
   Digit specification parser
   ========================= */

/* Parse strings like:
 *   "123", "1K", "10M", "2g", "132876K", "1e6", "3E7"
 *
 * Suffixes:
 *   K = 1_000
 *   M = 1_000_000
 *   G = 1_000_000_000
 *   T = 1_000_000_000_000
 *
 * Scientific notation: "<int>e<int>"
 *
 * On success: returns 0 and stores result in *out_digits (unsigned long).
 * On error:   returns non-zero and prints a message into stderr.
 */
static int parse_digit_spec(const char *spec, unsigned long *out_digits) {
    if (!spec) {
        fprintf(stderr, "Empty digits specification\n");
        return 1;
    }

    while (isspace((unsigned char)*spec)) spec++;  // trim leading spaces
    size_t len = strlen(spec);
    while (len > 0 && isspace((unsigned char)spec[len - 1])) len--; // trim trailing
    if (len == 0) {
        fprintf(stderr, "Empty digits specification\n");
        return 1;
    }

    /* Copy trimmed substring into a temp buffer */
    char buf[256];
    if (len >= sizeof(buf)) {
        fprintf(stderr, "Digits specification too long\n");
        return 1;
    }
    memcpy(buf, spec, len);
    buf[len] = '\0';

    /* 1) Scientific notation: "<int>e<int>" or "<int>E<int>" */
    char *e_ptr = strchr(buf, 'e');
    if (!e_ptr) {
        e_ptr = strchr(buf, 'E');
    }

    if (e_ptr) {
        *e_ptr = '\0';
        const char *mantissa_str = buf;
        const char *exp_str = e_ptr + 1;

        if (*mantissa_str == '\0' || *exp_str == '\0') {
            fprintf(stderr, "Invalid scientific notation: \"%s\"\n", spec);
            return 1;
        }

        unsigned long long mantissa = 0;
        unsigned long long exp = 0;

        if (sscanf(mantissa_str, "%llu", &mantissa) != 1) {
            fprintf(stderr, "Invalid mantissa in \"%s\"\n", spec);
            return 1;
        }
        if (sscanf(exp_str, "%llu", &exp) != 1) {
            fprintf(stderr, "Invalid exponent in \"%s\"\n", spec);
            return 1;
        }

        /* multiplier = 10^exp, with overflow checks */
        unsigned long long multiplier = 1;
        for (unsigned long long i = 0; i < exp; i++) {
            if (multiplier > ULLONG_MAX / 10ULL) {
                fprintf(stderr, "Exponent too large in \"%s\"\n", spec);
                return 1;
            }
            multiplier *= 10ULL;
        }

        if (mantissa > ULLONG_MAX / multiplier) {
            fprintf(stderr, "Digits value overflow for \"%s\"\n", spec);
            return 1;
        }
        unsigned long long value = mantissa * multiplier;

        if (value > ULONG_MAX) {
            fprintf(stderr,
                    "Too many digits (%llu), max supported is %lu\n",
                    value, (unsigned long)ULONG_MAX);
            return 1;
        }

        *out_digits = (unsigned long)value;
        return 0;
    }

    /* 2) Suffix-based notation: K, M, G, T */
    char last = buf[len - 1];
    unsigned long long multiplier = 1;
    if (last == 'k' || last == 'K' ||
        last == 'm' || last == 'M' ||
        last == 'g' || last == 'G' ||
        last == 't' || last == 'T') {

        switch (toupper((unsigned char)last)) {
            case 'K': multiplier = 1000ULL; break;
            case 'M': multiplier = 1000000ULL; break;
            case 'G': multiplier = 1000000000ULL; break;
            case 'T': multiplier = 1000000000000ULL; break;
        }
        buf[len - 1] = '\0';  // remove suffix
        len--;
        if (len == 0) {
            fprintf(stderr, "Missing number before suffix in \"%s\"\n", spec);
            return 1;
        }
    }

    unsigned long long base = 0;
    if (sscanf(buf, "%llu", &base) != 1) {
        fprintf(stderr, "Invalid number part in \"%s\"\n", spec);
        return 1;
    }

    if (base > ULLONG_MAX / multiplier) {
        fprintf(stderr, "Digits value overflow for \"%s\"\n", spec);
        return 1;
    }
    unsigned long long value = base * multiplier;

    if (value > ULONG_MAX) {
        fprintf(stderr,
                "Too many digits (%llu), max supported is %lu\n",
                value, (unsigned long)ULONG_MAX);
        return 1;
    }

    *out_digits = (unsigned long)value;
    return 0;
}

/* Get digits from command-line arguments.
 *
 * Supported forms:
 *   ./pi_chudnovsky               -> default (100000)
 *   ./pi_chudnovsky 12345
 *   ./pi_chudnovsky --calculate 1K
 *   ./pi_chudnovsky --digits 10M
 *   ./pi_chudnovsky -c 321
 *   ./pi_chudnovsky -d 132876K
 *   ./pi_chudnovsky -c 1e6
 */
static int get_digits_from_args(int argc, char **argv, unsigned long *out_digits) {
    const char *digit_spec = NULL;

    for (int i = 1; i < argc; i++) {
        const char *arg = argv[i];

        if (strcmp(arg, "--calculate") == 0 ||
            strcmp(arg, "-c") == 0 ||
            strcmp(arg, "--digits") == 0 ||
            strcmp(arg, "-d") == 0) {

            if (i + 1 >= argc) {
                fprintf(stderr, "Flag %s requires a value\n", arg);
                return 1;
            }
            digit_spec = argv[i + 1];
            i++; // skip value
        } else if (arg[0] != '-' && digit_spec == NULL) {
            /* First bare argument: treat as digits spec */
            digit_spec = arg;
        } else {
            /* Ignore other flags */
        }
    }

    /* Default if nothing given */
    if (digit_spec == NULL) {
        *out_digits = 100000UL;
        return 0;
    }

    return parse_digit_spec(digit_spec, out_digits);
}

/* =========================
   Chudnovsky binary split
   ========================= */

/*
 * Binary splitting for the Chudnovsky series.
 *
 * We compute P(a, b), Q(a, b), T(a, b) such that:
 *   π = (Q(0, N) * 426880 * sqrt(10005)) / T(0, N)
 *
 * a, b are term indices; P, Q, T are outputs (must be initialized by caller).
 */
static void binary_split(unsigned long a, unsigned long b,
                         mpz_t P, mpz_t Q, mpz_t T) {
    if (b - a == 1) {
        if (a == 0) {
            mpz_set_ui(P, 1);
            mpz_set_ui(Q, 1);
            mpz_set_ui(T, 13591409UL);
        } else {
            unsigned long k = a;

            /* P_k = (6k - 5)(2k - 1)(6k - 1) */
            mpz_t term1, term2, term3;
            mpz_init(term1);
            mpz_init(term2);
            mpz_init(term3);

            mpz_set_ui(term1, 6UL * k - 5UL);
            mpz_set_ui(term2, 2UL * k - 1UL);
            mpz_set_ui(term3, 6UL * k - 1UL);

            mpz_mul(P, term1, term2);
            mpz_mul(P, P, term3);

            /* Q_k = k^3 * C^3 / 24, where C = 640320
            * C^3 / 24 = 10939058860032000
            */
            mpz_t k3, c3_over_24;
            mpz_init(k3);
            mpz_init(c3_over_24);

            mpz_ui_pow_ui(k3, k, 3UL);  /* k^3 */

            /* Portable across 32- and 64-bit: load constant as a decimal string */
            mpz_set_str(c3_over_24, "10939058860032000", 10);

            mpz_mul(Q, k3, c3_over_24);

            /* T_k = (-1)^k * (13591409 + 545140134 k) * P_k */
            mpz_t val;
            mpz_init(val);
            mpz_set_ui(val, 545140134UL);
            mpz_mul_ui(val, val, k);
            mpz_add_ui(val, val, 13591409UL);

            if (k % 2 == 1) {
                mpz_neg(val, val);
            }

            mpz_mul(T, val, P);

            mpz_clear(term1);
            mpz_clear(term2);
            mpz_clear(term3);
            mpz_clear(k3);
            mpz_clear(c3_over_24);
            mpz_clear(val);
        }
    } else {
        unsigned long m = (a + b) / 2;

        mpz_t P1, Q1, T1;
        mpz_t P2, Q2, T2;
        mpz_init(P1); mpz_init(Q1); mpz_init(T1);
        mpz_init(P2); mpz_init(Q2); mpz_init(T2);

        binary_split(a, m, P1, Q1, T1);
        binary_split(m, b, P2, Q2, T2);

        /* T(a, b) = Q(m, b) * T(a, m) + P(a, m) * T(m, b) */
        mpz_t tmp;
        mpz_init(tmp);

        mpz_mul(T, Q2, T1);       /* T = Q2 * T1 */
        mpz_mul(tmp, P1, T2);     /* tmp = P1 * T2 */
        mpz_add(T, T, tmp);       /* T += tmp */

        /* Reuse P1 and Q1 as results: P = P1 * P2, Q = Q1 * Q2 */
        mpz_mul(P1, P1, P2);
        mpz_mul(Q1, Q1, Q2);

        mpz_set(P, P1);
        mpz_set(Q, Q1);

        mpz_clear(P1); mpz_clear(Q1); mpz_clear(T1);
        mpz_clear(P2); mpz_clear(Q2); mpz_clear(T2);
        mpz_clear(tmp);
    }
}

/* =========================
   Main
   ========================= */

int main(int argc, char **argv) {
    unsigned long digits;
    if (get_digits_from_args(argc, argv, &digits) != 0) {
        fprintf(stderr, "Usage examples:\n");
        fprintf(stderr, "  %s\n", argv[0]);
        fprintf(stderr, "  %s 12345\n", argv[0]);
        fprintf(stderr, "  %s --calculate 1K\n", argv[0]);
        fprintf(stderr, "  %s --digits 10M\n", argv[0]);
        fprintf(stderr, "  %s 1e6\n", argv[0]);
        return 1;
    }

    printf("Calculating pi to %lu digits (C + GMP + MPFR, Chudnovsky)...\n",
           digits);

    clock_t start = clock();

    /* Number of Chudnovsky terms, ~digits / 14 */
    unsigned long terms = digits / 14 + 1;

    mpz_t P, Q, T;
    mpz_init(P);
    mpz_init(Q);
    mpz_init(T);

    binary_split(0, terms, P, Q, T);

    /* Precision in bits: bits ~ digits * log2(10) + margin */
    double bits_per_digit = 3.321928094887362; /* log2(10) */
    double extra_bits = 256.0;
    mpfr_prec_t prec = (mpfr_prec_t)(digits * bits_per_digit + extra_bits);

    /* MPFR variables */
    mpfr_t sqrt10005, num, den, pi, scale, pi_scaled, pi_floor;
    mpfr_init2(sqrt10005, prec);
    mpfr_init2(num, prec);
    mpfr_init2(den, prec);
    mpfr_init2(pi, prec);
    mpfr_init2(scale, prec);
    mpfr_init2(pi_scaled, prec);
    mpfr_init2(pi_floor, prec);

    /* sqrt(10005) */
    mpfr_set_ui(sqrt10005, 10005UL, MPFR_RNDN);
    mpfr_sqrt(sqrt10005, sqrt10005, MPFR_RNDN);

    /* numerator = (Q * 426880) * sqrt(10005) */
    mpz_t Q_times_c;
    mpz_init(Q_times_c);
    mpz_mul_ui(Q_times_c, Q, 426880UL);
    mpfr_set_z(num, Q_times_c, MPFR_RNDN);
    mpfr_mul(num, num, sqrt10005, MPFR_RNDN);

    /* denominator = T */
    mpfr_set_z(den, T, MPFR_RNDN);

    /* pi = numerator / denominator */
    mpfr_div(pi, num, den, MPFR_RNDN);

    /* scale = 10^digits (as MPFR) */
    mpfr_ui_pow_ui(scale, 10UL, digits, MPFR_RNDN);

    /* pi_scaled = pi * 10^digits */
    mpfr_mul(pi_scaled, pi, scale, MPFR_RNDN);

    /* floor to truncate (no rounding) */
    mpfr_floor(pi_floor, pi_scaled);

    /* convert to integer */
    mpz_t pi_int;
    mpz_init(pi_int);
    mpfr_get_z(pi_int, pi_floor, MPFR_RNDN);

    clock_t end = clock();
    double elapsed = (double)(end - start) / (double)CLOCKS_PER_SEC;
    printf("Time: %.4f s\n", elapsed);

    /* Convert to string in base 10 */
    char *pi_str = mpz_get_str(NULL, 10, pi_int);
    size_t len = strlen(pi_str);

    /* We expect "3" + digits decimal digits. Ensure we have at least that. */
    size_t needed = digits + 1;
    if (len < needed) {
        /* Left-pad with zeros */
        size_t pad = needed - len;
        char *padded = (char *)malloc(needed + 1);
        if (!padded) {
            fprintf(stderr, "Out of memory while padding pi string\n");
            free(pi_str);
            goto cleanup;
        }
        memset(padded, '0', pad);
        memcpy(padded + pad, pi_str, len + 1); /* include '\0' */
        free(pi_str);
        pi_str = padded;
        len = needed;
    }

    /* Print as 3.<digits> */
    putchar(pi_str[0]);
    putchar('.');
    fwrite(pi_str + 1, 1, digits, stdout);
    putchar('\n');

cleanup:
    free(pi_str);

    mpz_clear(P);
    mpz_clear(Q);
    mpz_clear(T);
    mpz_clear(Q_times_c);
    mpz_clear(pi_int);

    mpfr_clear(sqrt10005);
    mpfr_clear(num);
    mpfr_clear(den);
    mpfr_clear(pi);
    mpfr_clear(scale);
    mpfr_clear(pi_scaled);
    mpfr_clear(pi_floor);

    return 0;
}
