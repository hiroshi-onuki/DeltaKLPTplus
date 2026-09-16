import argparse
import sys
import time
from parameters import Parameters
from sqisign import SQIsign
from ring_sqisign import RingSQIsign

TESTS = ("base", "ring")


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
    parser.add_argument(
        "--num-parties", "--num_parties",
        dest="num_parties",
        type=positive_integer,
        default=3,
        help="number of parties in the ring benchmark (default: 3)",
    )
    parser.add_argument(
        "--tests",
        nargs="+",
        choices=TESTS,
        default=TESTS,
        help="benchmarks to run (default: all)",
    )
    return parser.parse_args(argv)


def benchmark_base(inst, num_trials):
    t_keygen = 0
    t_sign = 0
    t_verify = 0
    for i in range(num_trials):
        # generate keys
        t0 = time.time()
        sk, pk = inst.Keygen()
        t_keygen += time.time() - t0

        # sign
        t0 = time.time()
        sign = inst.Sign(sk, pk, b"Test message")
        t_sign += time.time() - t0

        # verify
        t0 = time.time()
        assert inst.Verify(pk, b"Test message", sign)
        t_verify += time.time() - t0

        print(f"\r\033[2KBase signature trial {i+1}/{num_trials} completed.\r", end="")
    print(f"\r\033[2K", end="")
    return float(t_keygen / num_trials), float(t_sign / num_trials), float(t_verify / num_trials)

def benchmark_ring(inst, num_trials):
    t_keygen = 0
    t_sign = 0
    t_verify = 0
    for i in range(num_trials):
        # generate keys
        t0 = time.time()
        Pk, Sk = inst.Keygen()
        t_keygen += time.time() - t0

        # sign
        idx = randint(0, inst.n_parties - 1) # type: ignore
        t0 = time.time()
        sign = inst.Sign(Pk, Sk[idx], idx, b"Test message")
        t_sign += time.time() - t0

        # verify
        t0 = time.time()
        assert inst.Verify(Pk, b"Test message", sign)
        t_verify += time.time() - t0

        print(f"\r\033[2KRing signature trial {i+1}/{num_trials} completed.\r", end="")
    print(f"\r\033[2K", end="")
    return float(t_keygen / num_trials), float(t_sign / num_trials), float(t_verify / num_trials)

def main(argv=None):
    args = parse_args(argv)
    selected_tests = set(args.tests)
    tests = [test for test in TESTS if test in selected_tests]

    total_time = time.time()
    print(
        f"Running benchmarks ({', '.join(tests)}) with {args.num_trials} "
        "trials for each parameter set.\n"
    )
    for param in Parameters:
        f, c, lam = param["f"], param["c"], param["lam"]
        base_instance = None
        ring_instance = None
        if "base" in selected_tests:
            base_instance = SQIsign(f, c, lam)
        if "ring" in selected_tests:
            ring_instance = RingSQIsign(f, c, lam, args.num_parties)

        instance = base_instance if base_instance is not None else ring_instance
        parameter_text = f"Parameters: f={f}, c={c}, lam={lam}, e_rsp={instance.e_rsp}"
        if ring_instance is not None:
            parameter_text += f", num_parties={ring_instance.n_parties}"
        print(parameter_text)

        if "base" in selected_tests:
            avg_keygen, avg_sign, avg_verify = benchmark_base(base_instance, args.num_trials)
            print(
                f"Base SQIsign - Avg Keygen: {avg_keygen:.6f}s, "
                f"Avg Sign: {avg_sign:.6f}s, Avg Verify: {avg_verify:.6f}s"
            )
        if "ring" in selected_tests:
            avg_rkeygen, avg_rsign, avg_rverify = benchmark_ring(ring_instance, args.num_trials)
            print(
                f"Ring SQIsign - Avg Keygen: {avg_rkeygen:.6f}s, "
                f"Avg Sign: {avg_rsign:.6f}s, Avg Verify: {avg_rverify:.6f}s"
            )
            if "base" in selected_tests:
                print(
                    f"Ring SQIsign - Avg Keygen / n : {avg_rkeygen / ring_instance.n_parties:.6f}s, "
                    f"Avg Sign dummy : {(avg_rsign-avg_sign)/(ring_instance.n_parties-1):.6f}s, "
                    f"Avg Sing real : {avg_sign:.6f}s, "
                    f"Avg Verify / n : {avg_rverify / ring_instance.n_parties:.6f}s"
                )
        print()

    print(f"Total time for selected benchmarks: {time.time() - total_time:.2f}s")


# The modular Sage CLI may execute .sage files with __name__ == "sage.all".
if __name__ in ("__main__", "sage.all"):
    main()
