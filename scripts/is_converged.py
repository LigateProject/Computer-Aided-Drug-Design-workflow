from pathlib import Path as path
from argparse import ArgumentParser as cmd_line_parser
from enum import Enum as enum
from statistics import mean as mean
from csv import reader
from sys import stderr, exit


class status(enum):
    OK = "YES"
    KO = "NO"

    def __str__(self):
        return str(self.value)


def log_error(error: str):
    print(status.OK)  # to avoid an endless loop in the caller script
    print(f"Error: {error}", file=stderr)


if __name__ == "__main__":

    # we need to take the path to the edge path from the command line arguments
    parser = cmd_line_parser(
        prog="Convy",
        description="Check whether the simulation reached a stable value",
    )
    parser.add_argument(
        "datfile",
        type=path,
        default=path("."),
        help='Path to the target ".dat" file',
    )
    parser.add_argument(
        "--max_error",
        type=float,
        required=False,
        default=1.0,
        help="Maximum error tollerance for the measure",
    )
    parser.add_argument(
        "--time_limit",
        type=float,
        required=False,
        default=1000000.0,
        help="Maximum time that we want to simulate",
    )
    parser.add_argument(
        "--separator",
        type=str,
        required=False,
        default=" ",
        help="Field separator in the .dat file",
    )
    parser.add_argument(
        "--window_size",
        type=int,
        required=False,
        default=3,
        help="Number of row(s) in the consideration window",
    )
    args = parser.parse_args()

    # input validation
    if args.max_error < 0.0:
        log_error(f'The threshold "{args.max_error}" must be positive')
        raise ValueError("Negative error")
    if args.window_size < 1:
        log_error(f'The window size "{args.window_size}" must be strictly positive')
        raise ValueError("Negative window size")
    if args.time_limit < 0.0:
        log_error(f'The time limit "{args.time_limit}" must be positive!')
        raise ValueError("Negative error")

    # read the CSV table from the file
    x = []
    y = []
    e = []
    try:
        with open(args.datfile, "r") as infile:
            csv = reader(infile, delimiter=args.separator)
            for row in csv:
                r = [x for x in row if x]  # remove empty columns
                x.append(float(r[0]))
                y.append(float(r[1]))
                e.append(float(r[2]))
    except IOError as err:
        log_error(f'Unable to read from datfile "{args.datfile}"')
        raise err
    except IndexError as err:
        log_error(f'A row in "{args.datfile}" does not have enough columns')
        raise err
    except ValueError as err:
        log_error(f'A row in "{args.datfile}" have NaN as value')
        raise err

    # if we hit the time limit, we can just exit
    if x[-1] > args.time_limit:
        print(status.OK)
        exit(0)

    # if the standard deviation of the last line is above the threshold we can bail out
    last_error = e[-1]
    if last_error > args.max_error:
        print(status.KO)
        exit(0)

    # if we don't have enough data we can't say much
    if len(x) < args.window_size:
        print(status.KO)
        exit(0)

    # use a moving window to check for stability
    window = y[-args.window_size :]
    mu = mean(window)
    min_error = min(e[-args.window_size :])
    if any(abs(x - mu) > min_error for x in window):
        print(status.KO)
        exit(0)

    # if we reach this statement, the measure is stable
    print(status.OK)
