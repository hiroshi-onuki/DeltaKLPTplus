from parameters import Parameters
from sqisign import SQIsign
from quaternion import DeltaKLPT_plus

for param in Parameters:
    SQIsign_instance = SQIsign(param["e"], param["f"], param["lam"])
    print(SQIsign_instance.e_rsp)