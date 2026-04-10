# v0.3

- Add support for AVX512 instructions.
- Release AVX512 binary.
- Enable monomial orders: deglex, lexic, grevlexw
- Add support for ARM64 and NEON SIMD instructions.
- Release MacOS binary.

# v0.2
- Significant decrease in peak memory usage.
- Moderate speed-up in linear algebra phase for p > 2^16.
- Minor changes in CLI input parameters and log messages.
- Fixed bug when reducing certain large matrices for p < 2^8.
- Fixed bug in linear algebra phase when reducing a huge number of S-pairs.

# v0.1
- Support prime characteristics in the range (2^1, 2^31).
- Support grevlex and block-grevlex monomial orders.

