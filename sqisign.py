from sage.all import (
    is_prime,
    sqrt,
    kronecker,
    EllipticCurve,
)
import special_curve
import quaternion

class SQIsign:
    def __init__(self, p, e, f, lam):
        E0withEnd = special_curve.SpecialSuperSingularCurve(p, e, f)
        self.E0withEnd = E0withEnd
        Dmix = p**2 + 2
        while not is_prime(Dmix):
            Dmix += 2
        self.Dmix = Dmix
        self.lam = lam

    def Keygen(self):
        Isk, _ = quaternion.RandomFixedNormIdeal(self.E0withEnd.order, self.Dmix)
        Epk, Psk, Qsk = special_curve.IdealToIsogeny(self.E0withEnd, Isk)
        Ppk, Qpk = self._deterministic_torsion_basis(Epk, self.E0withEnd.e)

        sk = (Isk, Psk, Qsk)
        pk = Epk
        return sk, pk

    def Hash(self, msg):
        shake = SHAKE256.new(msg)
        return int.from_bytes(shake.read(self.lam // 8))

    def Sign(self, sk, pk, msg):
        Isk, Ppk, Qpk = sk
        Epk = pk
        c = self.Hash(msg)

        alpha = quaternion.SmallGenerator(I)
        alphaP, alphaQ = self.E0withEnd.quaternion_action(alpha) 
        assert (alphaP + c*alphaQ).is_zero(), "KernelToIdeal failed"
        return (alpha, R)


    @staticmethod
    def _deterministic_torsion_basis(E, e):
        """
        Compute a deterministic basis of E[2^e] for a supersingular Montgomery curve
        """
        F = E.base_ring()
        p = F.characteristic()
        i = sqrt(F(-1))
        A = E.a2()
        assert A != 0
        assert E == EllipticCurve(F, [0, A, 0, 1, 0]) # y^2 = x^3 + A*x^2 + x

        h = 0
        if A.is_square():
            # We need the denominator 1 + i*h to be a non-square in F_{p^2},
            # which holds iff its norm 1 + h^2 is a non-square in F_p.
            while True:
                h += 1
                xP = -1/(1 + i*h)*A
                if kronecker(1 + h**2, p) == 1:
                    continue
                if (xP**3 + A*xP**2 + xP).is_square():
                    break
        else:
            while True:
                h += 1
                xP = h*A
                if (xP**3 + A*xP**2 + xP).is_square():
                    break
        cofactor = (p + 1) // (2**e)
        P = cofactor * E.lift_x(xP)
        Q = cofactor * E.lift_x(-xP - A)
        assert P.weil_pairing(Q, 2**e)**(2**(e - 1)) == -1
        return P, Q
