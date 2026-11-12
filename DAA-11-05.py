"""Compare two OTU networks by the symmetric difference of their edges."""
import argparse
import time

from smg_core import build_network, changed_degrees
from smg_input import add_network_arguments, nonnegative_float, read_input_csv
from smg_output import threshold_context, write_changed_degrees, write_thresholds


def main(argv=None):
    parser = argparse.ArgumentParser(description="Compare networks using each group's adaptive threshold")
    parser.add_argument('-in1', required=True, help='First input CSV file path')
    parser.add_argument('-in2', required=True, help='Second input CSV file path')
    parser.add_argument('-o', '--output', required=True, help='Output CSV file path')
    parser.add_argument('-pos1', type=nonnegative_float, default=None,
                        help='Positive correlation threshold of file 1 (default: its adaptive threshold)')
    parser.add_argument('-pos2', type=nonnegative_float, default=None,
                        help='Positive correlation threshold of file 2 (default: its adaptive threshold)')
    add_network_arguments(parser)
    args = parser.parse_args(argv)

    networks = []
    try:
        for path, fixed in ((args.in1, args.pos1), (args.in2, args.pos2)):
            data = read_input_csv(path)
            networks.append(build_network(
                data.iloc[:, 0].to_numpy(), data.iloc[:, 1:].to_numpy(),
                angle=args.angle, multiplier=args.lam, fixed=fixed, column_normalize=args.column_normalize,
            ))
    except (ValueError, OSError) as exc:
        parser.error(str(exc))

    for i, network in enumerate(networks, 1):
        threshold = network.threshold
        print(f'Threshold file {i}: {threshold.value:.6f} ({threshold.kind}, {threshold_context(threshold)})')
    standard_file = write_thresholds(args.output, [network.threshold for network in networks])
    print(f'Standard OTUs and thresholds saved to {standard_file}')
    print('Processing direct comparison ...')
    names, degrees = changed_degrees(*networks)
    write_changed_degrees(args.output, names, degrees)
    print(f'Results successfully saved to {args.output}')


if __name__ == '__main__':
    start = time.perf_counter()
    main()
    print(f'run time: {time.perf_counter() - start:.2f} s')
