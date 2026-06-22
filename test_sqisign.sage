from sqisign import SQIsign
from quaternion import SmallGenerator

e = 248
f = 5
p = 2^e * f - 1  
lam = 128

SQIsign_instance = SQIsign(1)
E0withEnd = SQIsign_instance.E0withEnd

for _ in range(10):
    sk, pk = SQIsign_instance.Keygen()
    print("Key generation successful")
    
    sign = SQIsign_instance.Sign(sk, pk, b"Test message")
    print("Signing successful")

    assert SQIsign_instance.Verify(pk, b"Test message", sign)
    print("Verification successful")

    chl = randint(0, 2**SQIsign_instance.e_chl - 1)
    com, rsp = SQIsign_instance.Simulator(pk, chl)
    print("Simulation successful")


