import time
from parameters import Parameters
from sqisign import SQIsign
from ring_sqisign import RingSQIsign
from quaternion import RandomFixedNormIdeal, DeltaKLPT_plus

def average_iterations(inst, num_trials):
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

        _, _, iterations = DeltaKLPT_plus(Icom, IskIchl, 2, inst.e_rsp, inst.sec_lambda, count_iter=True)
        total_iterations += iterations
        print(f"\r\033[2KCount iteration trial {i+1}/{num_trials} completed. Iterations: {iterations}\r", end="")
    return float(total_iterations / num_trials)

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
    return float(t_keygen / num_trials), float(t_sign / num_trials), float(t_verify / num_trials)

num_trials = 2
num_parties = 3
total_time = time.time()
print(f"Running benchmarks with {num_trials} trials for each parameter set.\n")
for param in Parameters:
    e, f, lam = param["e"], param["f"], param["lam"]
    base_instance = SQIsign(e, f, lam)
    ring_instance = RingSQIsign(e, f, lam, num_parties)
    print(f"Parameters: e={e}, f={f}, lam={lam}, e_rsp={base_instance.e_rsp}, num_parties={ring_instance.n_parties}")
    avg_iter = average_iterations(base_instance, num_trials)
    avg_keygen_base, avg_sign_base, avg_verify_base = benchmark_base(base_instance, num_trials)
    avg_keygen_ring, avg_sign_ring, avg_verify_ring = benchmark_ring(ring_instance, num_trials)
    print(f"Average iterations for DeltaKLPT_plus: {avg_iter:.2f}")
    print(f"Base SQIsign - Avg Keygen: {avg_keygen_base:.6f}s, Avg Sign: {avg_sign_base:.6f}s, Avg Verify: {avg_verify_base:.6f}s")
    print(f"Ring SQIsign - Avg Keygen: {avg_keygen_ring:.6f}s, Avg Sign: {avg_sign_ring:.6f}s, Avg Verify: {avg_verify_ring:.6f}s\n")
print(f"Total time for all benchmarks: {time.time() - total_time:.2f}s")
