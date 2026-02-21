# MT1464 Assignment 1 - Global Hull Girder Mechanics

This repository contains scripts and outputs for MT1464 Assignment 1.

## Requirements

- `python3`
- `pandoc` (for `.docx` report generation)
- Python packages:
  - `matplotlib`
  - `pypdf`

Install Python packages with:

```bash
python3 -m pip install matplotlib pypdf
```

## Run The Solver

Default run:

```bash
python3 scripts/solve_assignment.py
```

Custom run:

```bash
python3 scripts/solve_assignment.py --studienummer 6470115 --output-dir outputs
```

## Generate The Verslag (.docx)

Default (uses `6470114` and `Jan Albert Driessen`):

```bash
python3 scripts/generate_docx_report.py
```

Custom:

```bash
python3 scripts/generate_docx_report.py --studienummer 6470115 --name "Jan Albert Driessen"
```

## Main Paths

- Source scripts: `scripts/`
- Assignment PDFs: `docs/pdfs/`
- Generated results: `outputs/`
- Generated report assets: `outputs/report_assets/`
