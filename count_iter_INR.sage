import argparse
import sys
import time
from parameters import Parameters
from sqisign import SQIsign
from ring_sqisign import RingSQIsign
from quaternion import RandomFixedNormIdeal, DeltaKLPT_plus

def positive_integer(value):
    value = int(value)
    if value <= 0:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return value


def parse_args(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    # Some Sage 10.8+ launchers require `--` and leave it in script arguments.
    if argv and argv[0] == "--":
        argv = argv[1:]

    parser = argparse.ArgumentParser(description="Benchmark SQIsign and Ring SQIsign")
    parser.add_argument(
        "--num-trials", "--num_trials",
        dest="num_trials",
        type=positive_integer,
        default=10,
        help="number of trials for each parameter set (default: 10)",
    )
    return parser.parse_args(argv)

def average_iterations(inst, num_trials, initial_reduce, shortest_alpha, search_prime, prime1mod4):
    total_iterations = 0
    for i in range(num_trials):

        # generate ideals
        Isk, _ = RandomFixedNormIdeal(inst.E0withEnd.order, inst.Dmix)
        while True:
            Ichl, _ = RandomFixedNormIdeal(inst.E0withEnd.order, 2**inst.sec_lambda)
            if Ichl.is_primitive():
                break
        IskIchl = Isk.intersection(Ichl)
        Icom, _ = RandomFixedNormIdeal(inst.E0withEnd.order, inst.Dmix)

        _, _, iterations = DeltaKLPT_plus(Icom, IskIchl, 2, inst.e_rsp, inst.sec_lambda, initial_reduce=initial_reduce, shortest_alpha=shortest_alpha, search_prime=search_prime, prime1mod4=prime1mod4, count_iter=True)
        total_iterations += iterations
        print(f"\r\033[2KCount iteration trial {i+1}/{num_trials} completed. Iterations: {iterations}\r", end="")
    print(f"\r\033[2K", end="")
    return float(total_iterations / num_trials)

def main(argv=None):
    args = parse_args(argv)
    total_time = time.time()
    for param in Parameters:
        f, c, lam = param["f"], param["c"], param["lam"]
        instance = SQIsign(f, c, lam)
        print(f"Parameters: f={f}, c={c}, lam={lam}, e_rsp={instance.e_rsp}")

        for options in [
            [False, False, False, False],
            [True, True, False, False],
            [True, True, True, False],
            [True, True, True, True]
        ]:
            initial_reduce, shortest_alpha, search_prime, prime1mod4 = options
            print(f"Options: initial_reduce={initial_reduce}, shortest_alpha={shortest_alpha}, search_prime={search_prime}, prime1mod4={prime1mod4}")
            avg_iter = average_iterations(instance, args.num_trials, initial_reduce, shortest_alpha, search_prime, prime1mod4)
            print(f"Average iterations for DeltaKLPT_plus: {avg_iter:.2f}")
        print()

    print(f"Total time for selected benchmarks: {time.time() - total_time:.2f}s")

# The modular Sage CLI may execute .sage files with __name__ == "sage.all".
if __name__ in ("__main__", "sage.all"):
    main()
