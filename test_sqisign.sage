from sqisign import SQIsign

SQIsign_instance = SQIsign(324, 3, 128)
E0withEnd = SQIsign_instance.E0withEnd

for _ in range(10):
    sk, pk = SQIsign_instance.Keygen()
    print("Key generation successful")
    
    sign = SQIsign_instance.Sign(sk, pk, b"Test message")
    print("Signing successful")

    assert SQIsign_instance.Verify(pk, b"Test message", sign)
    print("Verification successful")

    chl = randint(0, 2**SQIsign_instance.e_chl - 1) # type: ignore
    com, rsp = SQIsign_instance.Simulator(pk, chl)
    print("Simulation successful")


