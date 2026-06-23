from sqisign import SQIsign
import util

class RingSQIsign(SQIsign):
    def __init__(self, sec_level, n_parties):
        super().__init__(sec_level)
        self.n_parties = n_parties

    def Keygen(self):
        Pk = []
        Sk = []
        for _ in range(self.n_parties):
            sk, pk = super().Keygen()
            Pk.append(pk)
            Sk.append(sk)
        return Pk, Sk

    def Sign(self, Pk, sk, idx, message):
        n = self.n_parties
        Rsp = [None] * self.n_parties
        mPk = b''.join([util.j_invariant_to_bytes(pk) for pk in Pk])

        com_true, st_true = self.Commit() # the signer's commitment
        com = com_true
        chl_start = None

        # make (chl, com, rsp) for the parties except for the signer
        for i in range(1, n):
            chl = super().Hash(message + util.j_invariant_to_bytes(com) + mPk)
            if (idx + i) % n == 0:
                chl_start = chl
            com, rsp = super().Simulator(Pk[(idx + i) % n], chl)
            Rsp[(idx + i) % n] = rsp

        # make (chl, com, rsp) for the signer
        chl = super().Hash(message + util.j_invariant_to_bytes(com) + mPk)
        rsp = super().Respond(sk, st_true, chl)
        Rsp[idx] = rsp
        if idx == 0:
            chl_start = chl

        for rsp in Rsp:
            assert rsp is not None

        return chl_start, Rsp

    def Verify(self, Pk, message, sign):
        chl_first, Rsp = sign
        mPk = b''.join([util.j_invariant_to_bytes(pk) for pk in Pk])

        chl = chl_first
        for i in range(self.n_parties):
            com, is_cyclic = super().RecoverCommitment(Pk[i], chl, Rsp[i])
            chl = super().Hash(message + util.j_invariant_to_bytes(com) + mPk)
            if not is_cyclic:
                print("not cyclic in", i)
                return False
        return chl == chl_first