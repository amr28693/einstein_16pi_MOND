# einstein\_16pi\_MOND

Reproduction code and data for:

> Rodriguez, Anderson M. (2026). *The MOND Interpolation Function as a Bose–Einstein Partition Function over Gravitational Entropy Modes: A Complexified Extension of the Einstein–Hilbert Action.* Submitted to *Entropy* (MDPI). DOI: [pending]

---

## For one-click result

```bash
python3 mond_bose_einstein.py rar_galaxies.txt
```

This runs the full analysis and prints:

- Residual statistics (mean, std, RMS) vs. the 2,693 SPARC data points
- Occupation number table across the full θ range
- The crossover prediction: θ = ln 2, g_bar/a₀ = (ln 2)² = 0.4805

Expected output (key lines):

```
Parsed 2693 data points
Mean residual:   -0.0036 dex
Std residual:    0.1329 dex
RMS residual:    0.1329 dex
RAR intrinsic scatter (McGaugh+2016): 0.13 dex
--> CONSISTENT: residuals within known RAR scatter
crossover_theta: 0.6931 (= ln 2)
crossover_g_over_a0: 0.4805 (= (ln 2)²)
n_BE_at_crossover: 1.0
```

Results are also saved to `bose_einstein_results.json` and `bose_einstein_plotdata.csv`.

---

## If you want the figures

**Figure 1** (the headline RAR comparison, three panels):

```bash
python3 make_rar_figure.py
```

Requires `bose_einstein_plotdata.csv` — run `mond_bose_einstein.py` first if not already present.

Outputs: `mond_be_rar_figure.pdf` and `mond_be_rar_figure.png`

**Appendix figures** (projection mechanism visualizations):

```bash
python3 regen_appendix_figures.py
```

Outputs: `R_i_total_curvature.png`, `soliton_gradient_fields.png`, `GR16pi_projection_simulation.png`

---

## Dependencies

Python 3.8 or later. Install with:

```bash
pip install numpy pandas matplotlib scipy
```

No other dependencies or GPU required. This should run quickly on any modern laptop.

---

## Repository contents

| File | Description |
|------|-------------|
| `mond_bose_einstein.py` | Core analysis script. Loads data, computes BE prediction, outputs statistics. |
| `make_rar_figure.py` | Generates the three-panel RAR comparison figure (Figure 1 in paper). |
| `regen_appendix_figures.py` | Regenerates the three appendix projection figures. |
| `rar_galaxies.txt` | SPARC RAR dataset — 2,693 measurements across 153 galaxies (see citation below). |
| `bose_einstein_results.json` | Cached numerical results from the analysis. |
| `bose_einstein_plotdata.csv` | Per-point computed values (θ, n_BE, D_obs, D_pred, residuals). |

---

## Data

`rar_galaxies.txt` contains the radial acceleration relation dataset from:

> Lelli, F.; McGaugh, S.S.; Schombert, J.M.; Pawlowski, M.S. (2017). One Law to Rule Them All: The Radial Acceleration Relation of Galaxies. *Astrophysical Journal*, 836, 152. https://doi.org/10.3847/1538-4357/836/2/152

Original data publicly available at: http://astroweb.cwru.edu/SPARC/

Please cite Lelli et al. (2017) if you use this dataset.

---

## Citation

If you use this code or reproduce these results, please cite:

```
Rodriguez, A.M. (2026). The MOND Interpolation Function as a Bose–Einstein 
Partition Function over Gravitational Entropy Modes: A Complexified Extension 
of the Einstein–Hilbert Action. Entropy (MDPI). DOI: [pending]
```

---

## Demonstrated within this code:

The MOND interpolation function

    g_obs / g_bar = 1 / (1 - exp(-θ))     where θ = sqrt(g_bar / a₀)

is algebraically identical to a Bose–Einstein partition function:

    g_obs / g_bar = 1 + n_BE(θ)           where n_BE = 1 / (exp(θ) - 1)

The "dark matter" contribution at each radius is the mean occupation number
of gravitational entropy modes. At the crossover θ = ln 2, exactly one bit
of entropy excitation appears — the Landauer cost of the real-sector projection
in the complexified Einstein–Hilbert action.

The prediction has zero free parameters beyond a₀, which is inherited from
existing MOND phenomenology. Tested against 2,693 SPARC measurements across
153 galaxies, residual scatter = 0.133 dex, equal to the known intrinsic RAR
scatter (McGaugh et al. 2016).

---

## License

MIT License

Copyright (c) 2026 Anderson M. Rodriguez, (except SPARC and RAR data, as noted above)

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

Data (`rar_galaxies.txt`): publicly available from Lelli et al. (2017) /
SPARC database (http://astroweb.cwru.edu/SPARC/). Please cite accordingly.
