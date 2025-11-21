
# PI CALC

- Pi decimal calculators implemented in different languages. 

- Calculates π to arbitrary precision using the Chudnovsky algorithm. 

- Implementations included in C, C++, Rust and Python (plain and gmpy2-backed). 

- Default behaviour is to compute 100000 digits. 

- C/C++/Rust versions are high-performance (~100M digits on an i7 in ~2 minutes); pure Python will be very slow for large targets. 

- Use suffixes (K, M, G, T) or scientific notation (1e6, 1E6).

# Features
- Multiple implementations (C, C++, Rust, Python, Python + gmpy2)
- Human-friendly number parsing: suffixes like 1K, 10M, 5G and scientific notation
- Flags: --digits and --calculate accepted by the programs
- Outputs π to stdout (redirect as needed)

# Requirements
- C / C++: gcc / g++ and GMP + MPFR (and libgmpxx for C++)
- Rust: cargo + rust toolchain
- Python: Python 3.x; optional gmpy2 for the gmpy2-backed script
- Build tools (make or manual compile as shown below)

# Build & run — examples
- Note: all implementations accept a decimal count argument (default 100000), numeric suffixes (K/M/G/T) and scientific forms (1e6, 1E6). 
- They also accept --digits and --calculate flags (examples below use a small set of commands — the same patterns apply to the other scripts/binaries).

# C
Build:
  
 ```
 gcc -O3 pi_chudnovsky.c -o pi_chudnovsky -lgmp -lmpfr
 ```

Run:
 
```
./pi_chudnovsky 1K
```

# C++
Build:

```
g++ -O3 pi_chudnovsky.cpp -o pi_chudnovsky_cpp -lgmpxx -lgmp -lmpfr
```
Run:

```
./pi_chudnovsky_cpp 1K
```
# Rust
Build & run:

```
cargo run --release
```
Run with argument:

```
cargo run --release -- 1K
```

# Python (pure)
Run:

```
python3 pi_chudnovsky.py 1K
```

# Python (gmpy2-backed)
Install optional dependency:
 
```
python3 -m venv venv; source ./venv/bin/activate && \
pip install gmpy2
```
Run:

```
source ./venv/bin/activate && \
python3 pi_chudnovsky_gmpy2.py 1K
```

# Flags & argument formats
- Positional argument: number of digits to compute (defaults to 100000)
- --digits <N> or --calculate <N>
- Suffixes: K (thousand), M (million), G (billion), T (trillion) — case-insensitive
- Scientific notation: 1e6 or 1E6 accepted
- Very large values (G/T or multi-million+) will require a lot of RAM and time — use with caution.

# Performance & notes
- The C, C++ and Rust implementations are optimized (GMP/MPFR-backed) and are significantly faster than the pure-Python version.
- Tested on an i7 ~100M digits computed in ≈2 minutes for C/C++/Rust builds (your hardware and configuration will vary).
- Python pure implementation is slow for large digit counts. The gmpy2-backed Python script will be almost as fast as c/cpp/rust implementation but depends on your gmpy2 build and environment.
- Output is printed to stdout; redirect to a file for large outputs (e.g., ./pi_chudnovsky 1M > pi-1M.txt).

# Contributing
- You are free to improve/share/use any of the code
