# Frozen HMD data provenance

`data/ITAFemaleData.rda` is the frozen Human Mortality Database (HMD) Italy
female dataset used to generate the synthetic populations in the paper. It is
a documented `StMoMoData` object, analogous to `StMoMo::EWMaleData`, and can be
inspected with `?ITAFemaleData`. Its source snapshot was obtained on 10
February 2026. The object is intentionally bundled with the package so that a
clean pipeline execution performs no data downloads and requires no HMD
credentials.

The one-year schedules consumed by the simulator cover calendar years
1872--2022 (151 years). The generator re-indexes these to simulation periods
0--150. Any manuscript statement that the input ends in 2021 should therefore
be updated to 2022 to agree with the frozen input and the reported results.

## Archive validation

Before conversion to the package dataset, the source snapshot was checked
against the `ITA/` subtree of
`hmd_countries_20260210.zip` (SHA-256
`2277582e03c2aa5632bb8d7f5350bdb2a38216c54d02bdce69af3b1205bbd29b`).
All 92 files matched byte-for-byte, with zero missing, extra, or mismatched
files.

The two source files used to construct `ITAFemaleData` have these SHA-256
checksums:

- `STATS/Mx_1x1.txt`:
  `07ddeb64683310cc61d9686ce25b8a389a0d21621e25298e307380e6f1b5fbaa`
- `STATS/Exposures_1x1.txt`:
  `6c8469db3fcf3e0be9a5cb73c0270857909c1a1e68a5db5816fe3cd8f0e819f9`

The resulting compressed `data/ITAFemaleData.rda` file has SHA-256
`c3e2fdc7387d108354da4f1e9a94c98d98bdbf5861f26b3538705e12d39d3d7c`.

## Reconciliation with the paper run

The deterministic preparation recipe is retained in
`data-raw/ITAFemaleData.R`. It refuses to build the package dataset unless both
source-file checksums match the values above, and it performs no downloads.

Using seed 1, super-population exposure 94,500, and sub-population exposures
5,000 and 500, a fresh simulation from this snapshot is equal to the cached
paper `dt_0` object and has the same serialized SHA-256 digest:
`1ab8407ecc0513a81aad45829b424462116c74e3a225e0467a9911e1f277a3ce`.
