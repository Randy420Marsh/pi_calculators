#include <gmpxx.h>
#include <mpfr.h>

#include <iostream>
#include <string>
#include <cctype>
#include <climits>
#include <chrono>

/* =========================
   Small helpers
   ========================= */

static std::string trim(const std::string &s) {
    const char *ws = " \t\r\n";
    auto start = s.find_first_not_of(ws);
    if (start == std::string::npos) return "";
    auto end = s.find_last_not_of(ws);
    return s.substr(start, end - start + 1);
}

/* Parse a digit specification like:
 * "123", "1K", "10M", "2g", "132876K", "1e6", "3E7"
 *
 * Suffixes:
 *   K = 1_000 (10^3)
 *   M = 1_000_000 (10^6)
 *   G = 1_000_000_000 (10^9)
 *   T = 1_000_000_000_000 (10^12)
 *
 * Scientific notation: "<int>e<int>", e.g. "1e6".
 *
 * Returns true on success, false on error.
 * Result stored into out_digits.
 */
static bool parse_digit_spec(const std::string &spec, unsigned long &out_digits) {
    std::string s = trim(spec);
    if (s.empty()) {
        std::cerr << "Empty digits specification\n";
        return false;
    }

    // 1) Scientific notation: "<int>e<int>" or "<int>E<int>"
    auto pos_e = s.find_first_of("eE");
    if (pos_e != std::string::npos) {
        std::string mantissa_str = s.substr(0, pos_e);
        std::string exp_str      = s.substr(pos_e + 1);

        if (mantissa_str.empty() || exp_str.empty()) {
            std::cerr << "Invalid scientific notation: \"" << spec << "\"\n";
            return false;
        }

        unsigned long long mantissa = 0;
        unsigned long long exp = 0;
        try {
            mantissa = std::stoull(mantissa_str);
            exp      = std::stoull(exp_str);
        } catch (...) {
            std::cerr << "Invalid mantissa or exponent in \"" << spec << "\"\n";
            return false;
        }

        // multiplier = 10^exp with overflow checks
        unsigned long long multiplier = 1;
        for (unsigned long long i = 0; i < exp; ++i) {
            if (multiplier > ULLONG_MAX / 10ULL) {
                std::cerr << "Exponent too large in \"" << spec << "\"\n";
                return false;
            }
            multiplier *= 10ULL;
        }

        if (mantissa > ULLONG_MAX / multiplier) {
            std::cerr << "Digits value overflow for \"" << spec << "\"\n";
            return false;
        }
        unsigned long long value = mantissa * multiplier;

        if (value > ULONG_MAX) {
            std::cerr << "Too many digits (" << value
                      << "), max supported is " << (unsigned long)ULONG_MAX << "\n";
            return false;
        }

        out_digits = static_cast<unsigned long>(value);
        return true;
    }

    // 2) Suffix-based notation: K, M, G, T
    char last = s.back();
    unsigned long long multiplier = 1;

    if (last == 'k' || last == 'K' ||
        last == 'm' || last == 'M' ||
        last == 'g' || last == 'G' ||
        last == 't' || last == 'T') {

        switch (std::toupper(static_cast<unsigned char>(last))) {
            case 'K': multiplier = 1000ULL; break;
            case 'M': multiplier = 1000000ULL; break;
            case 'G': multiplier = 1000000000ULL; break;
            case 'T': multiplier = 1000000000000ULL; break;
        }
        s.pop_back(); // remove suffix
        s = trim(s);
        if (s.empty()) {
            std::cerr << "Missing number before suffix in \"" << spec << "\"\n";
            return false;
        }
    }

    unsigned long long base = 0;
    try {
        base = std::stoull(s);
    } catch (...) {
        std::cerr << "Invalid number part in \"" << spec << "\"\n";
        return false;
    }

    if (base > ULLONG_MAX / multiplier) {
        std::cerr << "Digits value overflow for \"" << spec << "\"\n";
        return false;
    }
    unsigned long long value = base * multiplier;

    if (value > ULONG_MAX) {
        std::cerr << "Too many digits (" << value
                  << "), max supported is " << (unsigned long)ULONG_MAX << "\n";
        return false;
    }

    out_digits = static_cast<unsigned long>(value);
    return true;
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
 *   ./pi_chudnovsky 1e6
 */
static bool get_digits_from_args(int argc, char **argv, unsigned long &out_digits) {
    std::string digit_spec;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg == "--calculate" || arg == "-c" ||
            arg == "--digits"    || arg == "-d") {

            if (i + 1 >= argc) {
                std::cerr << "Flag " << arg << " requires a value\n";
                return false;
            }
            digit_spec = argv[++i];
        } else if (arg.size() > 0 && arg[0] != '-' && digit_spec.empty()) {
            // First bare argument: treat as digits spec
            digit_spec = arg;
        } else {
            // ignore other flags
        }
    }

    if (digit_spec.empty()) {
        out_digits = 100000UL;  // default
        return true;
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
 */
static void binary_split(unsigned long a, unsigned long b,
                         mpz_class &P, mpz_class &Q, mpz_class &T) {
    if (b - a == 1) {
        if (a == 0) {
            P = 1;
            Q = 1;
            T = 13591409L;
        } else {
            unsigned long k = a;

            // P_k = (6k - 5)(2k - 1)(6k - 1)
            mpz_class term1(6UL * k - 5UL);
            mpz_class term2(2UL * k - 1UL);
            mpz_class term3(6UL * k - 1UL);
            P = term1 * term2 * term3;

            // Q_k = k^3 * C^3 / 24, where C = 640320
            // C^3 / 24 = 10939058860032000
            mpz_class k3 = mpz_class(k) * k * k;      // k^3
            mpz_class c3_over_24("10939058860032000"); // <-- fixed
            Q = k3 * c3_over_24;

            // T_k = (-1)^k * (13591409 + 545140134 k) * P_k
            mpz_class val = 545140134L;
            val *= k;
            val += 13591409L;
            if (k % 2 == 1) {
                val = -val;
            }

            T = val * P;
        }
    } else {
        unsigned long m = (a + b) / 2;

        mpz_class P1, Q1, T1;
        mpz_class P2, Q2, T2;

        binary_split(a, m, P1, Q1, T1);
        binary_split(m, b, P2, Q2, T2);

        // T(a, b) = Q(m, b) * T(a, m) + P(a, m) * T(m, b)
        T = Q2 * T1 + P1 * T2;

        // P(a, b) = P(a, m) * P(m, b)
        // Q(a, b) = Q(a, m) * Q(m, b)
        P = P1 * P2;
        Q = Q1 * Q2;
    }
}

/* =========================
   Main
   ========================= */

int main(int argc, char **argv) {
    unsigned long digits;
    if (!get_digits_from_args(argc, argv, digits)) {
        std::cerr << "Usage examples:\n"
                  << "  " << argv[0] << "\n"
                  << "  " << argv[0] << " 12345\n"
                  << "  " << argv[0] << " --calculate 1K\n"
                  << "  " << argv[0] << " --digits 10M\n"
                  << "  " << argv[0] << " 1e6\n";
        return 1;
    }

    std::cout << "Calculating pi to " << digits
              << " digits (C++ + GMP/MPFR, Chudnovsky)...\n";

    auto start = std::chrono::high_resolution_clock::now();

    // Chudnovsky terms (~14 digits per term)
    unsigned long terms = digits / 14 + 1;

    mpz_class P, Q, T;
    binary_split(0, terms, P, Q, T);

    // Precision in bits: bits ≈ digits * log2(10) + margin
    const double bits_per_digit = 3.321928094887362; // log2(10)
    const double extra_bits     = 256.0;
    mpfr_prec_t prec = static_cast<mpfr_prec_t>(digits * bits_per_digit + extra_bits);

    // MPFR variables
    mpfr_t sqrt10005, num, den, pi, scale, pi_scaled, pi_floor;
    mpfr_init2(sqrt10005, prec);
    mpfr_init2(num,       prec);
    mpfr_init2(den,       prec);
    mpfr_init2(pi,        prec);
    mpfr_init2(scale,     prec);
    mpfr_init2(pi_scaled, prec);
    mpfr_init2(pi_floor,  prec);

    // sqrt(10005)
    mpfr_set_ui(sqrt10005, 10005UL, MPFR_RNDN);
    mpfr_sqrt(sqrt10005, sqrt10005, MPFR_RNDN);

    // numerator = (Q * 426880) * sqrt(10005)
    mpz_class Q_times_c = Q * 426880UL;
    mpfr_set_z(num, Q_times_c.get_mpz_t(), MPFR_RNDN);
    mpfr_mul(num, num, sqrt10005, MPFR_RNDN);

    // denominator = T
    mpfr_set_z(den, T.get_mpz_t(), MPFR_RNDN);

    // pi = numerator / denominator
    mpfr_div(pi, num, den, MPFR_RNDN);

    // scale = 10^digits
    mpfr_ui_pow_ui(scale, 10UL, digits, MPFR_RNDN);

    // pi_scaled = pi * 10^digits
    mpfr_mul(pi_scaled, pi, scale, MPFR_RNDN);

    // floor to truncate (no rounding)
    mpfr_floor(pi_floor, pi_scaled);

    // convert to integer
    mpz_class pi_int;
    mpfr_get_z(pi_int.get_mpz_t(), pi_floor, MPFR_RNDN);

    auto end = std::chrono::high_resolution_clock::now();
    double elapsed =
        std::chrono::duration<double>(end - start).count();

    std::cout << "Time: " << elapsed << " s\n";

    // Convert to base-10 string
    std::string pi_str = pi_int.get_str(10);
    std::size_t len = pi_str.size();
    std::size_t needed = static_cast<std::size_t>(digits) + 1; // "3" + digits decimals

    if (len < needed) {
        // left-pad with zeros
        std::string padded(needed - len, '0');
        padded += pi_str;
        pi_str.swap(padded);
        len = needed;
    }

    // Print as 3.<digits>
    std::cout << pi_str[0] << '.';
    std::cout.write(pi_str.data() + 1, digits);
    std::cout << '\n';

    // Cleanup MPFR
    mpfr_clear(sqrt10005);
    mpfr_clear(num);
    mpfr_clear(den);
    mpfr_clear(pi);
    mpfr_clear(scale);
    mpfr_clear(pi_scaled);
    mpfr_clear(pi_floor);

    return 0;
}
