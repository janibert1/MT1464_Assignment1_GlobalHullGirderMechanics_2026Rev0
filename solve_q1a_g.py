"""Backwards-compatible entrypoint.
Prefer: python scripts/solve_assignment.py --studienummer 6470114 --output-dir outputs
"""
import argparse
from scripts.solve_assignment import main

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--studienummer', default='6470114')
    parser.add_argument('--output-dir', default='outputs')
    args = parser.parse_args()
    main(args.studienummer, args.output_dir)
