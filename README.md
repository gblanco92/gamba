# GamBa

`GamBa` is a fast program for computing Groebner bases in polynomial rings. 

`GamBa` features a high-performance C++ implementation of the Monte-Carlo F4 algorithm. It currently supports prime characteristics in the range $(2^1, 2^{31})$ together with the `grevlex` and the `block-grevlex` monomial orders.

## Benchmarks

`GamBa` is currently the **fastest** system for computing Groebner bases over *finite fields* on a *single thread*.

Here is a comparison with other implementations of the F4 algorithm. Systems are run with default parameters and the best of three timings is recorded. All times are in seconds.

> *Benchmark system*: AMD Ryzen 5 3600 @3.60Ghz 6-Core Processor (Zen 2), 64 GiB of dual-channel RAM @3200 MT/s. 
TurboCore = Off, Multithreading = Off, CPU_Mitigations = Off. (Ubuntu 22.04.4 LTS)

| $\ \ \ \ p = 32003$ | gamba v0.1  | magma v2.28  | msolve v0.6.8 | Grobner.jl v0.7.5 | FGb v1.68
|---|:---:|:---:|:---:|:---:|:---:|
| cyclic9  | **17.8** | 50.6 | 184.88 | 90.5 | 66.8 |
| cyclic10 | **852.7** | 3,363.2  | 29,102.8 | 10,775.0 | 7,233.2 |
| katsura14  | **59.4**  | 132.9  | 1,116.7  | 347.4 | 181.8 |
| katsura15  | **318.5**  | 924.2  | 9,207.15 | 2,195.3 | 1,255.5 |
| eco14  | **16.6** | 39.2 | 282.9 | 115.9 | 49.2 |
| eco15  | **81.0**  | 272.0 | 2,320.0 | 701.8 | 304.8 |
| noon8 | **3.9** | 14.0 | 48.2 | 52.3 | 15.8 |
| noon9 | **5.8** | 13.8 | 8.0 | 18.5 | 7.8 |
| CP(4,9,9) | **81.1** | 357.4 | 2,096.8 | 974.4 | 501.3 |
| CP(4,10,10) | **1,607.7** | 6,623.0 | 54,387.2 | $-$ | 10,573.2 |
| reimer8 | **7.7** | 21.0 | 26.4 | 22.9 | 17.0 |
| reimer9 | **315.0** | 981.9 | 1,333.0 | 755.3 | 726.3 |
| minrank(2,2,13,2) | **4,107.4** | 11,860.6 | $-$ | $-$ | 17,884.2 |
| minrank(3,2,9,2) | **1,404.8** | 6,640.4 | 48,338.8 | $-$ | 10,311.3 |
| bayes148 | **8.3** | 9.7 | 41.4 | 19.0 | 15.4 |
| game2 | **2.3** | 6.9 | 15.9 | 19.8 | 8.1 |
| jason210 | **2.3** | 5.9 | 3.2 | 3.9 | 3.0 |
| mayr42 | 19.3 | **4.2** | 51.4 | 24.1 | 27.5 |
| yang1 | 9.4 | **5.5** | 78.0 | 11.1 | 17.0 |

<br>

| $\ \ \ p = 2^{31}-1$  | gamba v0.1  | magma v2.28  | msolve v0.6.8 | Grobner.jl v0.7.5 | FGb v1.68
|---|:---:|:---:|:---:|:---:|:---:|
| cyclic9  | **26.1** | 56.8 | 53.2 | 90.5 | 65.6 |
| cyclic10 | **1,262.0** | 5,860.3 | 5,119.0  | 10,976.3 | 7,176.5 |
| katsura14  | **74.3** | 233.2  | 120.6  | 350.5 | 178.7 |
| katsura15  | **423.3** | 1,692.7  | 783.4 | 2,211.8 | 1,238.6 |
| eco14  | **19.6** | 88.7 | 50.0 | 118.1 | 48.1 |
| eco15  | **97.7** | 538.3 | 322.3 | 727.0 | 298.6 |
| noon8 | **5.4** | 18.8 | 11.0 | 53.2 | 15.6 |
| noon9 | 8.9 | 19.2 | **7.8** | 18.6 | 7.8 |
| CP(4,9,9) | **130.5** | 525.0 | 307.3 | 973.0 | 498.2 |
| CP(4,10,10) | **2,623.5** | 10,544.6 | 6,141.3 | $-$ | 10,555.7 |
| reimer8 | **9.2** | 17.4 | 16.3 | 22.9 | 16.4 |
| reimer9 | **380.8** | 935.6 | 745.9 | 766.2 | 692.5 |
| minrank(2,2,13,2) | **6,214.1** | 16,631.8 | 10,502.2 | $-$ | 17,858.7 |
| minrank(3,2,9,2) | **2461.0** | 9,915.0 | 5,979.1 | $-$ | 10,251.8 |

<sub>Entries with a dash mean that the computation either did not finish in less than 24 hours or it ran out of memory.</sub>

## Usage

To compute a Grobner basis with `GamBa` simply run

`./gamba -i input_file.txt -o output_file.txt`

See the [examples](examples/) folder for the self-explanatory input file format.

Run `./gamba --help` for a description of the other options.

## TODO

- Multithreading support
- Rational reconstruction
- Use AVX512 and ARM-Neon instruction sets
- Release MacOS binary
- Write library interface
- Change of order algorithms
- Implement Ore/Weyl algebras

## Requiriments

The `GamBa` binary requires a Linux x86-64 system with `glibc` version 2.24 or newer and a CPU supporting the AVX2 instruction set.

## Licensing

The public source code of `GamBa` is distributed under the GNU General Public License 3.0 (GPLv3), see the [LICENSE](LICENSE.md) document. Some components of `GamBa` will remain closed source for the time being.
