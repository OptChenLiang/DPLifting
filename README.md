# DPLifting

**DPLifting** (v1.4.1) is a standalone C/C++ library for **DPL/DPT hybrid coefficient lifting** on general knapsack sets. Optional capacity reduction (**+R**) is independent of the DPL/DPT backend. Prefer the C++ API `lifting()` or the C ABI `dplifting_lift_cover()`.

## Scope

- General upper bounds (bounded and unbounded variables)
- Sequential up/down lifting with **DPL** or **DPT** subproblem tables
- Three DPL/DPT controls: hard force, τ-threshold policy, or automatic feature map
- Optional **+R** (compile with `REDUCTION=1`; runtime ON/OFF/AUTO)
- No external solver dependency

**Maturity:** \(\sum w x \le b\) (`is_leq = 1`) is the production path.  
\(\sum w x \ge b\) (`is_leq = 0`) is supported and covered by tests.

## Requirements

- C++11 (GCC or Clang; MSVC untested)
- GNU Make

## Build

| Target | Command |
| ------ | ------- |
| Library | `make` → `libdplifting.so` |
| Tests | `make test` |
| Extended tests | `make test-all` |
| Example | `make example && ./example` |

```bash
make && make test
make clean && make REDUCTION=0 && make test REDUCTION=0   # build without +R
```

```bash
g++ -O2 my_app.cpp -Iinclude -L. -ldplifting -Wl,-rpath,'$ORIGIN' -lm -o my_app
make install   # default PREFIX=$HOME/.local
```

## Parameter model (v1.3)

Two **orthogonal** axes. Do not combine `MODE_AUTO` with `threshold` to mean “force DPL/DPT”.

### Axis 1 — DPL / DPT (`isdpl_mode`)

| Mode | Meaning |
| ---- | ------- |
| `DPLIFTING_MODE_DPL` / `_DPT` | **Manual hard force.** No mid-lift switch. |
| `DPLIFTING_MODE_THRESHOLD` | **Manual τ policy.** `threshold` = capacity \(\tau\). Initial and mid-lift use \(\tau\) vs \(\bar b^k=\min(b^k,U^k)\): \(\tau>\bar b^k\) → DPT, \(\tau<\bar b^k\) → DPL. |
| `DPLIFTING_MODE_AUTO` | **Default automatic.** Select from row features \((\rho_w,\beta,\bar u)\). **Ignores `threshold` for mode.** No mid-lift. |

Defaults for AUTO (override via `lift->rho_th` / `beta_th` / `u_bar_th` before the call; `0` → library defaults):

| Symbol | Default | Rule of thumb |
| ------ | ------- | ------------- |
| \(\rho_w=w_{\max}/w_{\min}\) | \(\rho_{\mathrm{th}}=6\) | \(\rho_w\ge\rho_{\mathrm{th}}\) → DPL |
| \(\beta=b/\mathrm{mean}(w)\) | \(\beta_{\mathrm{th}}=6\) | \(\beta<\beta_{\mathrm{th}}\) → DPL; else if narrow \(w\) → DPT |
| \(\bar u=\mathrm{mean}(u)\) | \(\bar u_{\mathrm{th}}=3\) | near-binary / small \(\bar u\) → DPL |
| Fuzzy band | 10% of thresholds | favors DPL |

Helpers: `dplifting_compute_features`, `dplifting_select_backend`, `dplifting_policy_default`.

### Axis 2 — Reduction +R (`lift->reduction_request`)

Requires compile-time `REDUCTION=1` (default in `make`). Always disabled if any variable is treated as unbounded.

Define the relative fill \(\lambda = b / \sum_i w_i u_i\) from the knapsack rhs \(b\) (`cap`) and all finite upper bounds.

| Value | Meaning |
| ----- | ------- |
| `DPLIFTING_RED_ON` / `_OFF` | **Manual** |
| `DPLIFTING_RED_AUTO` (0, default) | Enable iff \(\lambda \ge \lambda_{\mathrm{th}}\) **and** \(\bar b^0 > \tau\) |

Default \(\lambda_{\mathrm{th}}=0.1\) (`DPLIFTING_DEFAULT_LAMBDA_TH`). Override with `lift->lambda_th` before the call (`0` → library default).  
**If \(\lambda < 0.1\), +R AUTO turns reduction off** (tight rows; reduction rarely helps and adds overhead).

### Role of `threshold`

| Context | Role of `threshold` |
| ------- | ------------------- |
| `MODE_THRESHOLD` | \(\tau\) for initial + mid-lift |
| `RED_AUTO` | \(\tau\) for the \(\bar b^0 > \tau\) test (only after the \(\lambda\) gate) |
| `MODE_AUTO` / `_DPL` / `_DPT` | **Not** used to choose the backend |

If `threshold <= 0`, \(\tau=\beta_{\mathrm{th}}\cdot\mathrm{mean}(w)\).

### Quick recipes

```cpp
DPLifting lift = {};
// Default production: AUTO backend + AUTO +R
lifting(&lift, ..., /*threshold*/ 0.0, 0.0, DPLIFTING_MODE_AUTO);

// Hard DPT / DPL (benchmarks)
lifting(&lift, ..., 0.0, 0.0, DPLIFTING_MODE_DPT);
lifting(&lift, ..., 0.0, 0.0, DPLIFTING_MODE_DPL);

// Manual τ with mid-lift (e.g. τ = 200)
lifting(&lift, ..., 200.0, 0.0, DPLIFTING_MODE_THRESHOLD);

// Force +R off
lift.reduction_request = DPLIFTING_RED_OFF;
lifting(&lift, ..., 0.0, 0.0, DPLIFTING_MODE_AUTO);
```

## API

Header: `include/DPLifting.h` (shim: `dplifting_c.h`). Version: `DPLIFTING_VERSION` (`"1.4.1"`).

### `lifting` (C++)

```cpp
int lifting(
      DPLifting* lift,
      double* p, double* w, double* u, int* isuseub,
      double cap, int isSubCap,
      int* seed, int n_seed,
      int* liftingorder, int n_liftingorder,
      double* rhs,
      int isLeq, double* x, int n,
      double threshold, double duration, int isdpl_mode);
```

| Parameter | Meaning |
| --------- | ------- |
| `lift` | Workspace. Optionally set `reduction_request`, `rho_th`, `beta_th`, `u_bar_th`, `lambda_th` before the call. On success `lift->duration` is CPU seconds; `feat_lambda` holds \(\lambda=b/\sum w_i u_i\) (struct is cleared/rebuilt inside; policy fields are preserved). |
| `p` | Seed coefficients in → lifted coefficients out |
| `w`, `u` | Weights and upper bounds |
| `isuseub` | `1` = down-lift (fix at UB), `0` = up-lift |
| `cap` / `isSubCap` | Capacity \(b\), or residual if `isSubCap=1` |
| `seed` / `liftingorder` | Cover indices and remaining lift order |
| `rhs` | Lifted RHS (out) |
| `isLeq` | `1`: \(\le\) knapsack; `0`: \(\ge\) (experimental) |
| `x` | Optional fractional point; may be `NULL` |
| `threshold` | \(\tau\) for `MODE_THRESHOLD` and `RED_AUTO`; see table above |
| `duration` | Optional time limit (seconds); `<=0` = unlimited |
| `isdpl_mode` | `MODE_AUTO` / `MODE_THRESHOLD` / `MODE_DPL` / `MODE_DPT` |

**Return:** `1` success, `0` failure.

### `dplifting_lift_cover` (C)

Same semantics; returns `DPLIFTING_OK` (0) or `DPLIFTING_ERR_*`.

## Example

Knapsack set \(\mathcal{X}=\{x\in\mathbb{Z}_+^5: 8x_1+5x_2+4x_3+3x_4+5x_5\le 23,\ \ldots\}\).  
Seed \(2x_1+x_2\le 4\) on a restricted face with residual capacity \(18\); lift order \(\{3,4,5\}\).

```cpp
#include <DPLifting.h>
#include <cstdio>

int main() {
   const int n = 5;
   double p[] = {2.0, 1.0, 0.0, 0.0, 0.0};
   double w[] = {8.0, 5.0, 4.0, 3.0, 5.0};
   double u[] = {2.0, 3.0, 6.0, 5.0, 1.0};
   int isuseub[] = {0, 0, 0, 0, 1};
   int seed[] = {0, 1};
   int order[] = {2, 3, 4};
   double rhs = 0.0;

   DPLifting lift = {};
   if (!lifting(&lift, p, w, u, isuseub, 18.0, 1,
            seed, 2, order, 3, &rhs, 1, nullptr, n,
            /* threshold */ 0.0, /* duration */ 0.0, DPLIFTING_MODE_DPT))
      return 1;

   for (int i = 0; i < n; i++)
      if (p[i] > EPS_DPL)
         std::printf("%.4f*x_%d + ", p[i], i + 1);
   std::printf("<= %.4f  (%.4f s)\n", rhs, lift.duration);
   return 0;
}
```

Expected: coef `[2, 1, 1, 0.5, 1.5]`, rhs `5.5`.

```bash
make && make example && ./example
```

## Layout

```
DPLifting/
├── include/DPLifting.h
├── include/dplifting_c.h
├── src/DPLifting.cpp
├── examples/example.cpp
├── tests/test_dplifting.cpp
├── Makefile
├── CHANGELOG.md
└── LICENSE
```

## License

MIT — see [LICENSE](LICENSE). Copyright (c) 2026 Xintong Wang, Liang Chen, Yu-hong Dai.

## Citation

```
Xintong Wang et al. DPLifting: DPL/DPT hybrid lifting for general knapsack set.
https://github.com/OptChenLiang/DPLifting (version 1.4.1).
https://159.226.92.34:8000/wangxintong/dllifting
```
