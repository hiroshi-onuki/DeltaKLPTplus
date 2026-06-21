from sage.all import (
    ZZ,
    is_prime,
    sqrt,
    kronecker,
    EllipticCurve,
    vector,
    norm,
    inverse_mod,
)
import hashlib
import special_curve
import quaternion
import util

class SQIsign:
    def __init__(self, sec_level):

        # NIST security level 1 parameters (lam = 128)
        if sec_level == 1:
            # use Sage Integers (ZZ) so downstream exact arithmetic (e.g. the
            # pairing exponent in special_curve) is not turned into Python floats
            self.e = ZZ(248)
            self.f = ZZ(5)
            p = 2**self.e * self.f - 1
            assert is_prime(p)
            self.p = p
            self.E0withEnd = special_curve.SpecialSuperSingularCurve(p, self.e, self.f)
            Dmix = p**2 + 2
            while not is_prime(Dmix):
                Dmix += 2
            self.Dmix = Dmix    # the degree of phi_sk and phi_com, which satisfies the mixing property in the supersingular isogeny graph
            self.EISN_norm_bound = ZZ(2)**260
            self.e_chl = ZZ(128)
            self.e_rsp = ZZ(1050)
        else:
            raise ValueError("Unsupported security level")

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
        return c % 2**self.e_chl

    def Sign(self, sk, pk, msg):
        Isk, Msk = sk
        Epk = pk
        O0 = self.E0withEnd.order
        e = self.E0withEnd.e
        Ppk, Qpk = self._deterministic_torsion_basis(Epk, e)

        Icom, _ = quaternion.RandomFixedNormIdeal(self.E0withEnd.order, self.Dmix)
        Ecom, Pcom, Qcom = self.E0withEnd.IdealToIsogeny(Icom)

        chl = self.Hash(msg + util.field_element_to_bytes(Ecom.j_invariant(), self.e_chl//2) + util.field_element_to_bytes(Epk.j_invariant(), self.e_chl//2))
        a, b = vector([1, chl]) * Msk.inverse()
        Ichl = self.E0withEnd.KernelToIdeal(a, b, self.e_chl)
        IskIchl = Isk.intersection(Ichl)

        IcomIrsp, _ = quaternion.deltaKLPTforSign(Icom, IskIchl, 2, self.e_rsp, self.EISN_norm_bound) # tmp!
        N = norm(IcomIrsp) / 2**self.e_rsp
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
        e0 = self.e_chl + self.e_rsp - 4*e
        
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

        v0 = self.E0withEnd.IdealToKernel(Im0p2, e0)
        v0 = v0 * Msk % 2**e0
        c0 = v0[1] * inverse_mod(v0[0], 2**e0) % 2**e0
        assert c0 % 2**self.e_chl == chl

        print("The norm of I1m1 is", norm(Im1))
        Em1, Pm1, Qm1 = self.E0withEnd.IdealToIsogeny(Im1)
        Pm1d, Qm1d = self._deterministic_torsion_basis(Em1, e)
        Mm1 = util.BiDLP_matrix_power_two(Pm1, Qm1, Pm1d, Qm1d, e)
        v1dual = self.E0withEnd.IdealToKernel(Im1p2b, e)
        v1dual = v1dual * Mm1 % 2**e
        K1dual = v1dual[0] * Pm1d + v1dual[1] * Qm1d
        # for check
        K = 2**(e-e0) * (Ppk + c0 * Qpk)
        Echl = Epk.isogeny(K, model='montgomery', algorithm='factored').codomain()
        Echld = Em1.isogeny(K1dual, model='montgomery', algorithm='factored').codomain()
        assert Echld.j_invariant() == Echl.j_invariant()


        print("The norm of I2m2 is", norm(Im2))
        Em2, Pm2, Qm2 = self.E0withEnd.IdealToIsogeny(Im2)
        Pm2d, Qm2d = self._deterministic_torsion_basis(Em2, e)
        Mm2 = util.BiDLP_matrix_power_two(Pm2, Qm2, Pm2d, Qm2d, e)


    @staticmethod
    def _deterministic_torsion_basis(E, e):
        """
        Compute a deterministic basis of E[2^e] for a supersingular Montgomery curve
        """
        F = E.base_ring()
        p = F.characteristic()
        i = F.gen()
        assert i**2 == -1
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
