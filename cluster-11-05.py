"""Build one OTU network and export its connected components."""
import argparse
import time

from smg_core import build_network, connected_components
from smg_input import add_network_arguments, add_output_arguments, nonnegative_float, read_input_csv, validate_output
from smg_output import threshold_context, write_binary_clusters, write_clusters


def main(argv=None):
    parser = argparse.ArgumentParser(description="Cluster OTUs using the group's adaptive threshold")
    parser.add_argument('-i', '--input', required=True, help='Input CSV file')
    add_output_arguments(parser, clustering=True)
    parser.add_argument('-t', '--threshold', type=nonnegative_float,
                        help="Custom threshold (optional; default: the group's adaptive threshold)")
    add_network_arguments(parser)
    args = parser.parse_args(argv)
    try:
        validate_output(args.output, args.format)
        data = read_input_csv(args.input)
        network = build_network(
            data.iloc[:, 0].to_numpy(), data.iloc[:, 1:].to_numpy(), angle=args.angle,
            multiplier=args.lam, fixed=args.threshold, column_normalize=args.column_normalize,
            include_distances=args.format == 'excel' or args.distances,
        )
    except (ValueError, OSError) as exc:
        parser.error(str(exc))

    threshold = network.threshold
    print(f'Using {threshold.kind} threshold: {threshold.value:.6f} ({threshold_context(threshold)})')
    groups = connected_components(network.adjacency)
    try:
        if args.format == 'binary':
            write_binary_clusters(args.output, network, groups)
        else:
            write_clusters(args.output, network, groups)
    except (ValueError, OSError) as exc:
        parser.error(str(exc))
    print('\nClustering Results:')
    for i, members in enumerate(groups, 1):
        if len(members) >= 2:
            print(f'Cluster {i}: {len(members)} members')


if __name__ == '__main__':
    start = time.perf_counter()
    main()
    print(f'run time: {time.perf_counter() - start:.2f} s')
