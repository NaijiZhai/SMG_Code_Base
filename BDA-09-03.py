"""Compute sample-to-sample SMG distances on unscaled input abundances."""
import argparse
import time

from smg_core import distance_matrix
from smg_input import add_output_arguments, read_input_csv, validate_output
from smg_output import write_binary_distances, write_distances


def main(argv=None):
    parser = argparse.ArgumentParser(description='Sample-to-sample beta difference analysis')
    parser.add_argument('-i', '--input', required=True, help='Input CSV file')
    add_output_arguments(parser)
    args = parser.parse_args(argv)
    try:
        validate_output(args.output, args.format)
        data = read_input_csv(args.input)
        names = data.iloc[:, 0].to_numpy()
        distances = distance_matrix(data.iloc[:, 1:].to_numpy())
        if args.format == 'binary':
            write_binary_distances(args.output, names, distances)
        else:
            write_distances(args.output, names, distances)
    except (ValueError, OSError) as exc:
        parser.error(str(exc))
    print(f'Results saved to {args.output}')


if __name__ == '__main__':
    start = time.perf_counter()
    main()
    print(f'run time: {time.perf_counter() - start:.2f} s')
