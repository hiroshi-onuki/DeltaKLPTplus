from ring_sqisign import RingSQIsign

e = 248
f = 5
p = 2^e * f - 1  
lam = 128

SQIsign_instance = RingSQIsign(1, 5)

for _ in range(10):
    Pk, Sk = SQIsign_instance.Keygen()
    print("Key generation successful")
    
    idx = randint(0, SQIsign_instance.n_parties - 1) # type: ignore
    print(f"Signing with party index {idx}")

    sign = SQIsign_instance.Sign(Pk, Sk[idx], idx, b"Test message")
    print("Signing successful")

    assert SQIsign_instance.Verify(Pk, b"Test message", sign)
    print("Verification successful")


