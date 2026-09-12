from parameters import Parameters
from sqisign import SQIsign

n = 10

for param in Parameters:
    e, f, lam = param["e"], param["f"], param["lam"]
    SQIsign_instance = SQIsign(e, f, lam)
    print(f"Parameters: e={e}, f={f}, lam={lam}, e_rsp={SQIsign_instance.e_rsp}")
    for i in range(n):
        sk, pk = SQIsign_instance.Keygen()
        sign = SQIsign_instance.Sign(sk, pk, b"Test message")
        assert SQIsign_instance.Verify(pk, b"Test message", sign)
        print(f"\r\033[2K    {i+1}/{n} completed.\r", end="")
    print("All tests passed for this parameter set\n")


