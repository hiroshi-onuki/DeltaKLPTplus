from parameters import Parameters
from ring_sqisign import RingSQIsign

n = 10

for param in Parameters:
    e, f, lam = param["e"], param["f"], param["lam"]
    SQIsign_instance = RingSQIsign(e, f, lam, n_parties=3)
    print(f"Parameters: e={e}, f={f}, lam={lam}, e_rsp={SQIsign_instance.e_rsp}")
    for i in range(n):
        Pk, Sk = SQIsign_instance.Keygen()
        idx = randint(0, SQIsign_instance.n_parties - 1) # type: ignore
        sign = SQIsign_instance.Sign(Pk, Sk[idx], idx, b"Test message")
        assert SQIsign_instance.Verify(Pk, b"Test message", sign)
        print(f"\r\033[2K    {i+1}/{n} completed.\r", end="")
    print("All tests passed for this parameter set\n")


