# GamBa

### A fast program for computing Groebner bases in polynomial rings. ###

`GamBa` features a high-performance C++ implementation of the Monte-Carlo F4 algorithm. It supports prime characteristics in the range $(2^1, 2^{31})$ and multiple monomial orders.

> **Warning:** Some parts of `GamBa` will remain closed-source for the time being, so it is not currently possible to build `GamBa` from source. If you want to try out `GamBa`, use the provided binaries.

## Benchmarks

`GamBa` is currently the **fastest** system for computing Groebner bases over *finite fields* on a *single thread*.

Here is a comparison with other implementations of the F4 algorithm. Systems are run with default parameters and the best of three timings is recorded. All times are in seconds.

> *Benchmark system*: AMD EPYC 9175F @4.20 Ghz 16-Core Processor (Zen 5), 768 GiB DDR5 12 channel RAM @5200 MT/s.
TurboCore = On, Multithreading = Off, CPU_Mitigations = Off. (Debian 14)

| $p = 32003$ | gamba v0.3  | magma v2.28  | msolve v0.9.4 | axf4 v2b | FGb v1.68
|---|:---:|:---:|:---:|:---:|:---:|
| cyclic9  | **8.0** | 19.8 | 104.2 | 41.9 | 42.9 |
| cyclic10 | **389.6** | 1,284.7 | 15,373.8 | 4,147.8 | 4,851.4 |
| katsura14  | **28.5** | 52.8 | 625.4 | 102.8 | 122.4 |
| katsura15  | **162.9** | 355.7 | 4,875.0 | 661.7 | 830.5 |
| eco14  | **8.1** | 17.3 | 161.5 | 52.6 | 31.5 |
| eco15  | **41.8** | 120.4 | 1,304.5 | 399.3 | 198.1 |
| noon8 | **0.3** | 0.9 | 0.5 | 1.1 | 0.5 |
| noon9 | **2.8** | 6.4 | 4.3 | 11.8 | 4.3 |
| CP(4,9,9) | **39.7** | 135.5 | 1,168.4 | 407.1 | 308.8 |
| CP(4,10,10) | **791.3** | 2,466.1 | 28,780.2 | 8,539.3 | 6,850.9 |
| reimer8 | **4.1** | 8.9 | 16.0 | 1,944.2 | 9.8 |
| reimer9 | **152.4** | 386.0 | 692.0 | $-^*$ | 396.7 |
| minrank(2,2,13,2) | **2,046.5** | 4,334.3 | 81,869.6 | 3,253.3 | 11,702.8 |
| minrank(3,2,9,2) | **684.7** | 2401.2 | 25,894.8 | 8,797.6 | 6,935.6 |
| bayes148 | **4.3** | 5.1 | 26.2 | 67.1 | 8.5 |
| game2 | **1.1** | 2.6 | 6.4 | 5.4 | 4.3 |
| jason210 | **1.2** | 2.8 | 1.9 | 4.9 | 1.6 |
| mayr42 | 9.6 | **2.1** | 29.6 | 9.0 | 15.1 |
| yang1 | 4.6 | **2.2** | 47.1 | $-$ | 7.4 |

<br>

| $p = 2^{31}-1$  | gamba v0.3  | magma v2.28  | msolve v0.9.4 | axf4 v2b | FGb v1.68
|---|:---:|:---:|:---:|:---:|:---:|
| cyclic9  | **9.7** | 29.6 | 31.3 | 41.9 | 41.6 |
| cyclic10 | **471.6** | 2,975.1 | 3,369.4 | 4,153.3 | 4,929.7 |
| katsura14  | **33.0** | 123.6 | 76.5 | 102,8 | 120.1 |
| katsura15  | **191.4** | 869.9 | 519.3 | 664.9 | 831.1 |
| eco14  | **9.6** | 47.8 | 28.8 | 52.6 | 30.5 |
| eco15  | **49.6** | 275.8 | 191.2 | 398.6 | 195.2 |
| noon8 | **0.4** | 1.1 | 0.5 | 1.1 | 0.5 |
| noon9 | **3.3** | 9.3 | 3.7 | 11.7 | 4.1 |
| CP(4,9,9) | **47.3** | 276.7 | 195.7 | 407.3 | 304.4 |
| CP(4,10,10) | **947.1** | 5,497.6 | 4,113.8 | 8,567.0 | 6,818.1 |
| reimer8 | **4.7** | 9.5 | 10.2 | 1,971.1 | 9.3 |
| reimer9 | **174.4** | 414.9 | 340.8 | $-$ | 362.5 |
| minrank(2,2,13,2) | **2,419.3** | 8,875.2 | 7,053.7 | 3,736.7 | 11,645.0 |
| minrank(3,2,9,2) | **842.5** | 5,067.4 | 4,013.8 | 10,577.5 | 6,904.6 |

$^*$ <sub>Entries with a dash mean that the computation either did not finish in less than 24 hours or they ran out of memory.</sub>

## Usage

To compute a Grobner basis with `GamBa` simply run

`./gamba -i input_file.txt -o output_file.txt`

See the [examples](examples/) folder for the self-explanatory input file format.

Run `./gamba --help` for a description of the other options.

## TODO

- Multithreading support
- Rational reconstruction
- Write library interface
- Change of order algorithms
- Implement Ore/Weyl algebras
- Support computations over modules

## Requiriments

The `GamBa` binary requires a Linux x86-64 system with `glibc` version >=2.24 and a CPU supporting the AVX2/AVX512 instruction set. A `GamBa` binary for MacOS running on Apple Silicon is also available.

## Licensing

The public source code of `GamBa` is distributed under the GNU General Public License 3.0 (GPLv3), see the [LICENSE](LICENSE.md) document.
