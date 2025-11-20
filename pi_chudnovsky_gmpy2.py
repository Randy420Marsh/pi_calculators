#!/usr/bin/env python3
"""
High-precision π calculator using Chudnovsky + binary splitting.

- Python + gmpy2 (GMP/MPFR under the hood).
- Near C/Rust speed for large digit counts.
- Floor-truncates the last digit, matching your C/C++/Rust behavior.
"""

from __future__ import annotations

import sys
from typing import Tuple

try:
    import gmpy2
    from gmpy2 import mpz, mpfr, sqrt, get_context, floor, digits as mpz_digits
except ImportError as e:
    sys.stderr.write(
        "Error: gmpy2 is not installed.\n"
        "Install it with, for example:\n"
        "  pip install gmpy2\n"
    )
    raise SystemExit(1)


# Optional: relax Python's int->str limit (won't affect gmpy2.mpz, but harmless)
if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)


# =========================
# Digit specification parser
# =========================


def parse_digit_spec(spec: str) -> int:
    """
    Parse a digit specification like:
      "123", "1K", "10M", "2g", "132876K", "1e6", "3E7"

    Suffixes (case-insensitive):
      K = 1_000 (10^3)
      M = 1_000_000 (10^6)
      G = 1_000_000_000 (10^9)
      T = 1_000_000_000_000 (10^12)

    Scientific notation:
      "<int>e<int>", e.g. "1e6".

    Returns: number of digits as Python int (unbounded).
    Raises ValueError on invalid input.
    """
    s = spec.strip()
    if not s:
        raise ValueError("Empty digits specification")

    # 1) Scientific notation: "<int>e<int>" or "<int>E<int>"
    for idx, ch in enumerate(s):
        if ch in ("e", "E"):
            mantissa_str = s[:idx]
            exp_str = s[idx + 1 :]
            if not mantissa_str or not exp_str:
                raise ValueError(f"Invalid scientific notation: {spec!r}")
            mantissa = int(mantissa_str)
            exp = int(exp_str)
            if exp < 0:
                raise ValueError(f"Negative exponent not supported in {spec!r}")
            value = mantissa * (10 ** exp)
            if value <= 0:
                raise ValueError(f"Digits must be positive: {spec!r}")
            return value

    # 2) Suffix-based notation: K, M, G, T
    last = s[-1]
    multiplier = 1
    if last in "kKmMgGtT":
        if last in "kK":
            multiplier = 1_000
        elif last in "mM":
            multiplier = 1_000_000
        elif last in "gG":
            multiplier = 1_000_000_000
        elif last in "tT":
            multiplier = 1_000_000_000_000
        s = s[:-1].strip()
        if not s:
            raise ValueError(f"Missing number before suffix in {spec!r}")

    base = int(s)
    if base <= 0:
        raise ValueError(f"Digits must be positive: {spec!r}")

    return base * multiplier


def get_digits_from_args(argv: list[str]) -> int:
    """
    Read digits from CLI arguments.

    Supported forms:
      python pi_chudnovsky_gmpy2.py            -> default (100000)
      python pi_chudnovsky_gmpy2.py 12345
      python pi_chudnovsky_gmpy2.py --calculate 1K
      python pi_chudnovsky_gmpy2.py --digits 10M
      python pi_chudnovsky_gmpy2.py -c 321
      python pi_chudnovsky_gmpy2.py -d 132876K
      python pi_chudnovsky_gmpy2.py 1e6
    """
    digit_spec: str | None = None
    args = argv[1:]  # skip program name

    i = 0
    while i < len(args):
        arg = args[i]
        if arg in ("--calculate", "-c", "--digits", "-d"):
            if i + 1 >= len(args):
                raise ValueError(f"Flag {arg!r} requires a value")
            digit_spec = args[i + 1]
            i += 2
        elif not arg.startswith("-") and digit_spec is None:
            # First bare argument treated as digits spec
            digit_spec = arg
            i += 1
        else:
            # Ignore other flags for now
            i += 1

    # Default if nothing given
    if digit_spec is None:
        return 100_000

    return parse_digit_spec(digit_spec)


# =========================
# Chudnovsky binary split
# =========================


def binary_split(a: int, b: int) -> Tuple[mpz, mpz, mpz]:
    """
    Binary splitting for the Chudnovsky series.

    Compute P(a, b), Q(a, b), T(a, b) such that:
      π = (Q(0, N) * 426880 * sqrt(10005)) / T(0, N)
    """
    if b - a == 1:
        if a == 0:
            P = mpz(1)
            Q = mpz(1)
            T = mpz(13591409)
            return P, Q, T
        else:
            k = mpz(a)

            # P_k = (6k - 5)(2k - 1)(6k - 1)
            term1 = mpz(6) * k - 5
            term2 = mpz(2) * k - 1
            term3 = mpz(6) * k - 1
            P = term1 * term2 * term3

            # Q_k = k^3 * C^3 / 24, where C = 640320
            # C^3 / 24 = 10939058860032000
            k3 = k * k * k
            c3_over_24 = mpz("10939058860032000")
            Q = k3 * c3_over_24

            # T_k = (-1)^k * (13591409 + 545140134 k) * P_k
            val = mpz(545140134) * k + 13591409
            if int(k) % 2 == 1:
                val = -val
            T = val * P

            return P, Q, T
    else:
        m = (a + b) // 2
        P1, Q1, T1 = binary_split(a, m)
        P2, Q2, T2 = binary_split(m, b)

        # P(a, b) = P(a, m) * P(m, b)
        P = P1 * P2

        # Q(a, b) = Q(a, m) * Q(m, b)
        Q = Q1 * Q2

        # T(a, b) = Q(m, b) * T(a, m) + P(a, m) * T(m, b)
        T = Q2 * T1 + P1 * T2

        return P, Q, T


# =========================
# Main computation
# =========================


def compute_pi_digits(digits: int) -> str:
    """
    Compute π to `digits` decimal places as a string "3.<digits>".

    Uses:
      - Chudnovsky + binary splitting in gmpy2.mpz
      - mpfr sqrt and division
      - extra precision and truncation (floor) to stabilize last digit
    """
    if digits <= 0:
        raise ValueError("digits must be positive")

    # Number of Chudnovsky terms (~14 digits per term)
    terms = digits // 14 + 1

    _P, Q, T = binary_split(0, terms)

    # Precision in bits: bits ≈ digits * log2(10) + safety margin
    bits_per_digit = 3.321928094887362  # log2(10)
    extra_bits = 256.0
    prec_bits = int(digits * bits_per_digit + extra_bits)

    ctx = get_context()
    ctx.precision = prec_bits

    # π = (Q * 426880 * sqrt(10005)) / T
    sqrt_10005 = sqrt(mpfr(10005))
    numerator = mpfr(Q) * 426880 * sqrt_10005
    denominator = mpfr(T)
    pi = numerator / denominator

    # Scale by 10^digits and truncate (floor) instead of round.
    scale = mpfr(10) ** digits
    pi_scaled = pi * scale

    # floor -> mpz, truncation (pi > 0 so floor == trunc)
    pi_floor = floor(pi_scaled)
    pi_int = mpz(pi_floor)

    # Convert to decimal string in base 10
    pi_str = mpz_digits(pi_int, 10)  # avoids Python int->str limits

    # Ensure at least digits + 1 characters ("3" + digits decimals)
    if len(pi_str) < digits + 1:
        pi_str = pi_str.rjust(digits + 1, "0")

    int_part = pi_str[0]
    frac_part = pi_str[1 : 1 + digits]
    return f"{int_part}.{frac_part}"


def main(argv: list[str]) -> int:
    try:
        digits = get_digits_from_args(argv)
    except ValueError as e:
        prog = argv[0] if argv else "pi_chudnovsky_gmpy2.py"
        sys.stderr.write(f"Error: {e}\n")
        sys.stderr.write("Usage examples:\n")
        sys.stderr.write(f"  {prog}\n")
        sys.stderr.write(f"  {prog} 12345\n")
        sys.stderr.write(f"  {prog} --calculate 1K\n")
        sys.stderr.write(f"  {prog} --digits 10M\n")
        sys.stderr.write(f"  {prog} 1e6\n")
        return 1

    print(f"Calculating π to {digits} digits (Python + gmpy2, Chudnovsky)...")

    import time

    start = time.perf_counter()
    pi_str = compute_pi_digits(digits)
    elapsed = time.perf_counter() - start

    print(f"Time: {elapsed:.6f} s")
    print(pi_str)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))

