# Package validation

Validation performed on 2026-08-18:

- MATLAB Code Analyzer parsed every `.m` file in the package without a syntax error.
- The MOSAIC solver SHA-256 is `8c3202da9b17ad6020f264203e0d872f8d84ac1dc07f19dede4bffce9efeb1c3` in all three subpackages.
- The PGL2021 callable solver hashes match across all three subpackages.
- The ZW2024 solver SHA-256 is `12e9d2ffb483f9266740643bd46c062b67f54f75d2161b94332141ed35d80bd6` in all three subpackages.
- The embedded truth MAT SHA-256 is `436c82138d4d40e4f6c90c455f96ddb13a43a2a776379e4886c2f5b6406e8aca` in all three subpackages.
- The selected p=5000 candidate identifiers are MOSAIC `(407,420,421)`, PGL2021 `(269,294,319)`, and ZW2024 `(593,482,592)`.
- The embedded candidate thresholds match the p=5000 component-F1 validation selections used by the conference result lineage.
- MATLAB sources contain no machine-specific user paths.

Run `VERIFY_PACKAGE.m` for the static table, path, and cross-subpackage solver checks. Run each subpackage's `RUN_QUICK_TEST.m` for a short end-to-end execution check in the target MATLAB installation.
