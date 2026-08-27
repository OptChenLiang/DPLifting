# Changelog

## [1.4.1] - 2026-08-27

### Changed

- **+R AUTO:** disable reduction when \(\lambda = b/\sum_i w_i u_i < \lambda_{\mathrm{th}}\)
  (default \(\lambda_{\mathrm{th}}=0.1\)). Otherwise keep \(\bar b^0 > \tau\) gate.
- New fields: `lift->lambda_th`, `lift->feat_lambda`; `DPLiftingFeatures.lambda`,
  `DPLiftingPolicy.lambda_th`; macro `DPLIFTING_DEFAULT_LAMBDA_TH`.

## [1.4.0] - 2026-08-20

### Changed

- **Rename project** `DLLifting` → `DPLifting`; algorithm labels **DL→DPL**, **DP→DPT**.
- Public macros: `DPLIFTING_MODE_DPL` / `_DPT`; parameter `isdpl_mode`; library `libdplifting.so`.

## [1.3.0] - 2026-08-12

### Changed

- **Orthogonal controls** for hybrid lifting (documented in README):
  - **DPL/DPT:** `MODE_DPL` / `MODE_DPT` (hard force); `MODE_THRESHOLD` (manual \(\tau\) vs
    \(\bar b^k\), mid-lift); `MODE_AUTO` (feature map on \(\rho_w,\beta,\bar u\); default).
  - **+R:** `reduction_request` = `RED_ON` / `RED_OFF` / `RED_AUTO` (\(\bar b^0>\tau\)).
- `threshold` is \(\tau\) only for `MODE_THRESHOLD` and `RED_AUTO` — not for `MODE_AUTO`.
- Helpers: `dplifting_policy_default`, `dplifting_compute_features`,
  `dplifting_select_backend`, `Lifting_bar_b`; Multiply mid-lift unified via `Lifting_add_chunk`.

## [1.2.2] - 2026-08-06

### Changed

- `duration` argument to `lifting()` is now an optional **time limit** (seconds):
  `<= 0` means unlimited; if the limit is exceeded during sequential lifting, `lifting()` returns 0.
  Measured runtime is still written to `lift->duration`.
- Cleared library `-Wall` unused-parameter/variable warnings.

## [1.2.1] - 2026-08-06

### Changed

- `make test` always runs geq (`>=`) suite; its failures affect the process exit code.
- Library: replace `assert` failure paths with `return 0` error codes in Up/Down/Lifting;
  propagate failure from `lifting()`; C ABI returns `DPLIFTING_ERR_INTERNAL`.
- Release builds keep `-DNDEBUG`; debug `printf` remains behind `DPLIFTING_DEBUG`.
- README: full markdown parameter table for `lifting(...)`.

## [1.2.0] - 2026-08-06

### Changed

- Public-release hygiene: library no longer calls `exit()` on realloc failure; debug
  `printf` is gated behind `DPLIFTING_DEBUG`; release builds use `-DNDEBUG`.
- `make test` runs the core (`<=`) suite only and fails the process on any core failure.
- `make test-all` runs extended `>=` and mixed suites (`./test_dplifting --all`).
- README documents `threshold` / `isdpl_mode`, full `lifting` / `dplifting_lift_cover`
  signatures, and marks `>=` knapsacks as experimental.
- `make install` installs both `DPLifting.h` and `dplifting_c.h`; default `PREFIX=$(HOME)/.local`.
- Single translation unit: C ABI lives in `src/DPLifting.cpp`.
- Version macros: `DPLIFTING_VERSION` (`1.2.0`).

### Fixed

- Removed dead realloc/`exit(0)` path in `Lifting_Mergesort`.

## [1.1.0] - 2026-05-28

### Added

- `isdpl_mode` on `lifting()` / `dplifting_lift_cover()` (`AUTO` / `DPL` / `DPT`).

## [1.0.0] - 2026-05

### Added

- Initial standalone release.
