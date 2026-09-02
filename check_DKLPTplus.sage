from parameters import Parameters
from sqisign import SQIsign
from quaternion import RandomFixedNormIdeal, DeltaKLPT_plus


def average_iterations(SQIsign_instance, num_trials):
    Isk, _ = RandomFixedNormIdeal(SQIsign_instance.E0withEnd.order, SQIsign_instance.Dmix)
    while True:
        Ichl, _ = RandomFixedNormIdeal(SQIsign_instance.E0withEnd.order, 2**lam)
        if Ichl.is_primitive():
            break
    IskIchl = Isk.intersection(Ichl)
    Icom, _ = RandomFixedNormIdeal(SQIsign_instance.E0withEnd.order, SQIsign_instance.Dmix)

    total_iterations = 0
    for i in range(num_trials):
        _, _, iterations = DeltaKLPT_plus(Icom, IskIchl, 2, SQIsign_instance.e_rsp, SQIsign_instance.sec_lambda, count_iter=True)
        total_iterations += iterations
        print(f"Trial {i+1}/{num_trials} completed. Iterations: {iterations}")
    return total_iterations / num_trials

for param in Parameters:
    e, f, lam = param["e"], param["f"], param["lam"]
    SQIsign_instance = SQIsign(e, f, lam)
    print(f"Parameters: e={e}, f={f}, lam={lam}")
    print(SQIsign_instance.e_rsp)
    print(f"Average iterations for DeltaKLPT_plus: {average_iterations(SQIsign_instance, 10)}")
