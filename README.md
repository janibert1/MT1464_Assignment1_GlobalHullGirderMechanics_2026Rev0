# MT1464 Assignment 1 – Global Hull Girder Mechanics

This repository contains the assignment PDFs and a reproducible Python solver for:
- **Question 1a–1g** (with a dedicated re-check report for 1a–1f), and
- **Question 2a–2f**, and
- **Question 3a–3d**.

The current default setup is for studienummer **6470114**.

## Repository structure

- `scripts/solve_assignment.py` → main solver implementation.
- `solve_q1a_g.py` → backward-compatible entrypoint.
- `outputs/` → latest generated answer sheets, data, and plots.
- `legacy_outputs/` → earlier generated output files retained for traceability.
- `*.pdf` → assignment and reference documents.

## Run

```bash
python scripts/solve_assignment.py --studienummer 6470114 --output-dir outputs
```

Generated files:
- `outputs/answer_q1a_g.md`
- `outputs/check_q1a_f.md`
- `outputs/answer_q2a_f.md`
- `outputs/answer_q3a_d.md`
- `outputs/q3a_plots.svg`
- `outputs/q1a_g_data.csv`
- `outputs/q1a_g_plots.svg`
- `outputs/q1f_compare.svg`

## Notes

- The script derives digits `a..g` from the supplied studienummer.
- Q1 uses digit `g` according to Table 2/3 mappings.
- Q2 uses digit `f` for `t_p,eq` according to Table 4.

- Q1 bending moment now includes a crane-induced point moment jump at the crane location: `M_crane = SWL * reach = 1.5 * 30 = 45 MNm` (positive jump).
