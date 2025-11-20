use rug::{Float, Integer};
use rug::ops::Pow;
use std::time::Instant;

/// Accepts arbitrary digit counts: 1, 123, 7615236, etc.
///
/// Accepts suffixes (case-insensitive):
///   K = 1_000
///   M = 1_000_000
///   G = 1_000_000_000
///   T = 1_000_000_000_000
///
/// Also accepts scientific notation: 1e6, 3E7, etc.
///
/// Works with:
///   cargo run --release                  (defaults to 100_000 digits)
///   cargo run --release -- 12345
///   cargo run --release -- --calculate 1K
///   cargo run --release -- --digits 10M

/// Binary splitting for the Chudnovsky series
///
/// We compute P(a, b), Q(a, b), T(a, b) such that:
///   π = (Q(0, N) * 426880 * sqrt(10005)) / T(0, N)
fn binary_split(a: u64, b: u64) -> (Integer, Integer, Integer) {
    if b - a == 1 {
        // Base case: single term with index a
        if a == 0 {
            let p = Integer::from(1);
            let q = Integer::from(1);
            let t = Integer::from(13591409);
            (p, q, t)
        } else {
            let k = a;

            // P_k = (6k - 5)(2k - 1)(6k - 1)
            let term1 = Integer::from(6 * k - 5);
            let term2 = Integer::from(2 * k - 1);
            let term3 = Integer::from(6 * k - 1);
            let p = term1 * term2 * term3;

            // Q_k = k^3 * C^3 / 24, where C = 640320
            // C^3 / 24 = 10939058860032000
            let c3_over_24 = Integer::from(10939058860032000u64);
            let q = Integer::from(k).pow(3u32) * c3_over_24;

            // T_k = (-1)^k * (13591409 + 545140134 k) * P_k
            let val = Integer::from(13591409u64 + 545140134u64 * k);
            let t = if k % 2 == 1 {
                // use -val (moves val) to avoid NegIncomplete
                -val * &p
            } else {
                val * &p
            };

            (p, q, t)
        }
    } else {
        // Recursive split
        let m = (a + b) / 2;

        // Left half
        let (mut p1, mut q1, t1) = binary_split(a, m);
        // Right half
        let (p2, q2, t2) = binary_split(m, b);

        // Merge:
        // P(a, b) = P(a, m) * P(m, b)
        // Q(a, b) = Q(a, m) * Q(m, b)
        // T(a, b) = Q(m, b) * T(a, m) + P(a, m) * T(m, b)

        // Compute T into a fresh Integer
        let mut t = Integer::from(&q2 * &t1);
        t += Integer::from(&p1 * &t2);

        // Reuse p1 and q1 as the result to reduce allocations
        p1 *= p2; // p1 = p1 * p2
        q1 *= q2; // q1 = q1 * q2

        (p1, q1, t)
    }
}

/// Parse a digit specification like:
/// "123", "1K", "10M", "2g", "132876K", "1e6", "3E7"
///
/// Supports suffixes:
///   K = 1_000 (10^3)
///   M = 1_000_000 (10^6)
///   G = 1_000_000_000 (10^9)
///   T = 1_000_000_000_000 (10^12)
///
/// Also supports scientific notation: "<int>e<int>", e.g. "1e6".
///
/// Returns number of digits as u32 (max ~4.29e9).
fn parse_digit_spec(spec: &str) -> Result<u32, String> {
    let s = spec.trim();
    if s.is_empty() {
        return Err("Empty digits specification".to_string());
    }

    // 1) Scientific notation: "<int>e<int>" or "<int>E<int>"
    if let Some(pos) = s.find(|c| c == 'e' || c == 'E') {
        let (mantissa_str, exp_str_with_e) = s.split_at(pos);
        let exp_str = &exp_str_with_e[1..]; // skip 'e' or 'E'

        if mantissa_str.is_empty() || exp_str.is_empty() {
            return Err(format!("Invalid scientific notation: \"{}\"", spec));
        }

        let mantissa: u64 = mantissa_str
            .parse()
            .map_err(|_| format!("Invalid mantissa in \"{}\"", spec))?;
        let exp: u32 = exp_str
            .parse()
            .map_err(|_| format!("Invalid exponent in \"{}\"", spec))?;

        let multiplier = 10u64
            .checked_pow(exp)
            .ok_or_else(|| format!("Exponent too large in \"{}\"", spec))?;

        let value = mantissa
            .checked_mul(multiplier)
            .ok_or_else(|| format!("Digits value overflow for \"{}\"", spec))?;

        if value > u32::MAX as u64 {
            return Err(format!(
                "Too many digits ({}), max supported is {}",
                value,
                u32::MAX
            ));
        }

        return Ok(value as u32);
    }

    // 2) Suffix-based notation: K, M, G, T
    let s_upper = s.to_uppercase();
    let (num_str, multiplier): (&str, u64) = if s_upper.ends_with('K') {
        (&s_upper[..s_upper.len() - 1], 1_000)
    } else if s_upper.ends_with('M') {
        (&s_upper[..s_upper.len() - 1], 1_000_000)
    } else if s_upper.ends_with('G') {
        (&s_upper[..s_upper.len() - 1], 1_000_000_000)
    } else if s_upper.ends_with('T') {
        (&s_upper[..s_upper.len() - 1], 1_000_000_000_000)
    } else {
        (s_upper.as_str(), 1)
    };

    if num_str.is_empty() {
        return Err(format!("Missing number before suffix in \"{}\"", spec));
    }

    let base: u64 = num_str
        .parse()
        .map_err(|_| format!("Invalid number part in \"{}\"", spec))?;

    let value: u64 = base
        .checked_mul(multiplier)
        .ok_or_else(|| format!("Digits value overflow for \"{}\"", spec))?;

    if value > u32::MAX as u64 {
        return Err(format!(
            "Too many digits ({}), max supported is {}",
            value,
            u32::MAX
        ));
    }

    Ok(value as u32)
}

/// Get digits from CLI arguments.
/// Supported forms:
///   - <prog>                 -> default (100000)
///   - <prog> 12345
///   - <prog> --calculate 1K
///   - <prog> --digits 10M
///   - <prog> -c 321
///   - <prog> -d 132876K
///   - <prog> -c 1e6
///   - <prog> -d 1E6
fn get_digits_from_args() -> Result<u32, String> {
    let mut args = std::env::args().skip(1); // skip program name
    let mut digit_spec: Option<String> = None;

    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--calculate" | "-c" | "--digits" | "-d" => {
                let value = args
                    .next()
                    .ok_or_else(|| format!("Flag {} requires a value", arg))?;
                digit_spec = Some(value);
            }
            // First bare argument that is not a flag: treat as digits spec
            _ if !arg.starts_with('-') && digit_spec.is_none() => {
                digit_spec = Some(arg);
            }
            // Ignore any other flags for now
            _ => {}
        }
    }

    // Default if nothing given
    let default_digits: u32 = 100_000;

    match digit_spec {
        Some(spec) => parse_digit_spec(&spec),
        None => Ok(default_digits),
    }
}

fn main() {
    // Read digits from CLI
    let digits_u32 = match get_digits_from_args() {
        Ok(d) => d,
        Err(e) => {
            eprintln!("Error: {}", e);
            eprintln!("Usage examples:");
            eprintln!("  cargo run --release");
            eprintln!("  cargo run --release -- 12345");
            eprintln!("  cargo run --release -- --calculate 1K");
            eprintln!("  cargo run --release -- --digits 10M");
            eprintln!("  cargo run --release -- 1e6");
            std::process::exit(1);
        }
    };

    let digits: usize = digits_u32 as usize;

    println!(
        "Calculating π to {} digits (Rust + Rug, Chudnovsky)...",
        digits
    );
    let start = Instant::now();

    // Each Chudnovsky term yields ~14 digits
    let terms = digits / 14 + 1;
    let (_p, q, t) = binary_split(0, terms as u64);

    // Compute required precision in bits:
    // bits ≈ digits * log2(10) + safety_margin
    let bits_per_digit = 3.321_928_094_887_362_f64; // log2(10)
    let extra_bits = 256.0; // safety margin
    let precision: u32 = (digits as f64 * bits_per_digit + extra_bits) as u32;

    // π = (Q * 426880 * sqrt(10005)) / T
    let sqrt_10005 = Float::with_val(precision, 10005).sqrt();

    // Numerator as Integer first, then convert once to Float
    let q_times_c = Integer::from(426880) * &q;
    let numerator = Float::with_val(precision, q_times_c) * sqrt_10005;
    let denominator = Float::with_val(precision, &t);
    let pi = numerator / denominator;

    // Scale by 10^digits using Integer -> Float (explicit, safe)
    let scale_int = Integer::from(10).pow(digits_u32);
    let scale = Float::with_val(precision, &scale_int);
    let pi_scaled = &pi * scale;

    // Truncate instead of round: floor first, then convert to Integer
    let pi_floor = pi_scaled.floor();
    let pi_integer = pi_floor
        .to_integer()
        .expect("Failed to convert floor(pi * 10^digits) to Integer");

    println!("Time: {:.4}s", start.elapsed().as_secs_f64());

    // Convert to string and insert decimal point after first digit: 3.<digits>
    let mut pi_str = pi_integer.to_string();

    // Ensure at least digits + 1 characters ("3" + digits decimals)
    if pi_str.len() <= digits {
        pi_str = format!("{:0>width$}", pi_str, width = digits + 1);
    }

    let (int_part, frac_part) = pi_str.split_at(1);
    println!("{}.{}", int_part, &frac_part[..digits]);
}

