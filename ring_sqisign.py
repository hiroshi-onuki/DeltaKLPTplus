from sqisign import SQISign

class RingSQISign(SQISign):
    def __init__(self, sec_level, n_parties):
        super().__init__(sec_level)
        Pk = []
        Sk = []
        for i in range(n_parties):
            sk, pk = self.Keygen()
            Pk.append(pk)
            Sk.append(sk)
        self.Pk = Pk
        self.Sk = Sk
        self.n_parties = n_parties
