from sage.all import (
    ZZ,
    is_prime,
    sqrt,
    kronecker,
    EllipticCurve,
    vector,
    norm,
)
import hashlib
import special_curve
import quaternion
import util

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
        Epk, Psk, Qsk = self.E0withEnd.IdealToIsogeny(Isk)
        Ppk, Qpk = self._deterministic_torsion_basis(Epk, self.E0withEnd.e)
        Msk = util.BiDLP_matrix_power_two(Psk, Qsk, Ppk, Qpk, self.E0withEnd.e)
        sk = (Isk, Msk)
        pk = Epk
        return sk, pk

    def Hash(self, msg):
        h = hashlib.sha256(msg)
        c = util.bytes_to_integer(h.digest())
        return c % 2**self.lam

    def Sign(self, sk, pk, msg):
        Isk, Msk = sk
        Epk = pk
        O0 = self.E0withEnd.order
        e = self.E0withEnd.e
        Ppk, Qpk = self._deterministic_torsion_basis(Epk, e)

        Icom, _ = quaternion.RandomFixedNormIdeal(self.E0withEnd.order, self.Dmix)
        Ecom, Pcom, Qcom = self.E0withEnd.IdealToIsogeny(Icom)

        c = self.Hash(msg + util.field_element_to_bytes(Ecom.j_invariant(), self.lam//2) + util.field_element_to_bytes(Epk.j_invariant(), self.lam//2))
        a, b = Msk.transpose() * vector([1, c])
        Ichl = self.E0withEnd.KernelToIdeal(a, b, self.lam)
        IskIchl = Isk.intersection(Ichl)

        print("Starting deltaKLPTforSign...")
        IcomIrsp, _ = quaternion.deltaKLPTforSign(Icom, IskIchl, 2, 1050, 2**260) # tmp!
        print("deltaKLPTforSign completed.")
        N = norm(IcomIrsp) / 2**1050
        Icom_d = IcomIrsp + O0 * N
        assert Icom.right_order().isomorphism_to(Icom_d.right_order()) != None
        O = IcomIrsp.right_order()
        Od = IskIchl.right_order()
        alpha = O.isomorphism_to(Od, conjugator=True)
        Iall = IskIchl * alpha.inverse() * IcomIrsp.conjugate() * alpha
        n = norm(Iall) / alpha.reduced_norm()
        alpha *= ZZ(sqrt(n))     # scale alpha so that norm(Iall) = norm(alpha)
        assert Iall == O0 * alpha
        assert alpha/2 not in O0
        e0 = self.lam + 1050 - 4*e
        
        Im0p2 = Iall + O0 * 2**e0
        
        Im1Im1p2f = Iall + O0 * 2**(e0 + 2*e) * norm(Isk)
        Im1 = Im1Im1p2f + O0 * 2**(e0 + e) * norm(Isk)
        Im1, beta, nIm1 = quaternion.EquivalentRandomPrimeIdeal(Im1)
        Im1p2f = Im1Im1p2f * (beta.conjugate() / nIm1) + O0 * 2**e
        Im1p2b = O0 * beta.conjugate() + O0 * 2**e

        Im2Im2p2f = Iall + O0 * 2**(e0 + 4*e) * norm(Isk)
        Im2 = Im2Im2p2f + O0 * 2**(e0 + 3*e) * norm(Isk)
        Im2, beta, nIm2 = quaternion.EquivalentRandomPrimeIdeal(Im2)
        Im2p2f = Im2Im2p2f * (beta.conjugate() / nIm2) + O0 * 2**e
        Im2p2b = O0 * beta.conjugate() + O0 * 2**e

        assert Isk.intersection(Im0p2).right_order().isomorphism_to(Im1.intersection(Im1p2b).right_order()) != None
        assert Im1Im1p2f.right_order().isomorphism_to(Im2.intersection(Im2p2b).right_order()) != None
        assert Im2Im2p2f.right_order().isomorphism_to(Icom.right_order()) != None


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
