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
import utilities.discrete_log

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
        Epk, (Psk, Qsk) = self._normalize_curve(Epk, (Psk, Qsk))
        Ppk, Qpk = self._deterministic_torsion_basis(Epk, self.E0withEnd.e)
        Msk = util.BiDLP_matrix_power_two(Psk, Qsk, Ppk, Qpk, self.E0withEnd.e)
        sk = (Isk, Msk)
        pk = Epk
        return sk, pk

    def Commit(self, pk, msg):
        Icom, _ = quaternion.RandomFixedNormIdeal(self.E0withEnd.order, self.Dmix)
        Ecom, _, _ = self.E0withEnd.IdealToIsogeny(Icom)
        com = Ecom
        st = Icom
        return com, st
    
    def Respond(self, sk, st, chl):
        Isk, Msk = sk
        Icom = st
        O0 = self.E0withEnd.order
        e = self.E0withEnd.e

        # make the ideal corresponding to chl
        a, b = vector([1, chl]) * Msk.inverse()
        Ichl = self.E0withEnd.KernelToIdeal(a, b, self.e_chl)
        IskIchl = Isk.intersection(Ichl)

        IcomIrsp, _ = quaternion.deltaKLPTforSign(Icom, IskIchl, 2, self.e_rsp, self.EISN_norm_bound)
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
        Im1, beta, _ = quaternion.EquivalentRandomPrimeIdeal(Im1)
        Im1p2f = Im1Im1p2f * (beta.conjugate() / (2**(e0 + e)*norm(Isk))) + O0 * 2**e
        Im1p2b = O0 * beta.conjugate() + O0 * 2**e

        Im2Im2p2f = Iall + O0 * 2**(e0 + 4*e) * norm(Isk)
        Im2 = Im2Im2p2f + O0 * 2**(e0 + 3*e) * norm(Isk)
        Im2, beta, _ = quaternion.EquivalentRandomPrimeIdeal(Im2)
        Im2p2f = Im2Im2p2f * (beta.conjugate() / (2**(e0 + 3*e)*norm(Isk))) + O0 * 2**e
        Im2p2b = O0 * beta.conjugate() + O0 * 2**e

        v0 = self.E0withEnd.IdealToKernel(Im0p2, e0)
        v0 = v0 * Msk % 2**e0
        c0 = v0[1] * inverse_mod(v0[0], 2**e0) % 2**e0
        assert c0 % 2**self.e_chl == chl
        c0without_chl = (c0 - chl) // 2**self.e_chl

        Em1, Pm1, Qm1 = self.E0withEnd.IdealToIsogeny(Im1)
        Em1, (Pm1, Qm1) = self._normalize_curve(Em1, (Pm1, Qm1))
        Pm1d, Qm1d = self._deterministic_torsion_basis(Em1, e)
        Mm1 = util.BiDLP_matrix_power_two(Pm1, Qm1, Pm1d, Qm1d, e)
        v1dual = self.E0withEnd.IdealToKernel(Im1p2b, e)
        v1dual = v1dual * Mm1 % 2**e
        K1dual = v1dual[0] * Pm1d + v1dual[1] * Qm1d
        if v1dual[0] % 2 == 0:
            evalP = Pm1d
        else:
            evalP = Qm1d
        phi = Em1.isogeny(K1dual, model='montgomery', algorithm='factored')
        K1 = phi(evalP)
        assert K1.order() == 2**e
        Em1b = phi.codomain()
        Em1b, (K1,) = self._normalize_curve(Em1b, (K1,))
        Pm1b, Qm1b = self._deterministic_torsion_basis(Em1b, e)
        a, b = utilities.discrete_log.BiDLP_power_two(K1, Pm1b, Qm1b, e, None)
        if a % 2 == 0:
            c1b = a * inverse_mod(b, 2**e) % 2**e
            isP1b = True
        else:
            c1b = b * inverse_mod(a, 2**e) % 2**e
            isP1b = False
        v1f = self.E0withEnd.IdealToKernel(Im1p2f, e)
        v1f = v1f * Mm1 % 2**e
        a, b = v1f
        if a % 2 == 0:
            c1f = a * inverse_mod(b, 2**e) % 2**e
            isP1f = True
        else:
            c1f = b * inverse_mod(a, 2**e) % 2**e
            isP1f = False

        Em2, Pm2, Qm2 = self.E0withEnd.IdealToIsogeny(Im2)
        Em2, (Pm2, Qm2) = self._normalize_curve(Em2, (Pm2, Qm2))
        Pm2d, Qm2d = self._deterministic_torsion_basis(Em2, e)
        Mm2 = util.BiDLP_matrix_power_two(Pm2, Qm2, Pm2d, Qm2d, e)
        v2dual = self.E0withEnd.IdealToKernel(Im2p2b, e)
        v2dual = v2dual * Mm2 % 2**e
        K2dual = v2dual[0] * Pm2d + v2dual[1] * Qm2d
        if v2dual[0] % 2 == 0:
            evalP = Pm2d
        else:
            evalP = Qm2d
        phi = Em2.isogeny(K2dual, model='montgomery', algorithm='factored')
        K2 = phi(evalP)
        assert K2.order() == 2**e
        Em2b = phi.codomain()
        Em2b, (K2,) = self._normalize_curve(Em2b, (K2,))
        Pm2b, Qm2b = self._deterministic_torsion_basis(Em2b, e)
        a, b = utilities.discrete_log.BiDLP_power_two(K2, Pm2b, Qm2b, e, None)
        if a % 2 == 0:
            c2b = a * inverse_mod(b, 2**e) % 2**e
            isP2b = True
        else:
            c2b = b * inverse_mod(a, 2**e) % 2**e
            isP2b = False
        v2f = self.E0withEnd.IdealToKernel(Im2p2f, e)
        v2f = v2f * Mm2 % 2**e
        a, b = v2f
        if a % 2 == 0:
            c2f = a * inverse_mod(b, 2**e) % 2**e
            isP2f = True
        else:
            c2f = b * inverse_mod(a, 2**e) % 2**e
            isP2f = False

        return (c0without_chl, c1b, c1f, c2b, c2f, isP1b, isP1f, isP2b, isP2f)

    def Hash(self, msg):
        h = hashlib.sha256(msg)
        c = util.bytes_to_integer(h.digest())
        return c % 2**self.e_chl

    def Sign(self, sk, pk, msg):
        com, st = self.Commit(pk, msg)
        chl = self.Hash(msg + util.j_invariant_to_bytes(com) + util.j_invariant_to_bytes(pk))
        rsp = self.Respond(sk, st, chl)
        c0without_chl = rsp[0]
        c0 = c0without_chl * 2**self.e_chl + chl
        rsp = (c0,) + rsp[1:]
        return rsp

    def Verify(self, pk, msg, sign):
        e = self.E0withEnd.e
        e0 = self.e_chl + self.e_rsp - 4*e
        c0, c1b, c1f, c2b, c2f, isP1b, isP1f, isP2b, isP2f = sign
        Epk = pk
        Ppk, Qpk = self._deterministic_torsion_basis(Epk, self.E0withEnd.e)
        K = 2**(e - e0) * (Ppk + c0 * Qpk)
        E = Epk.isogeny(K, model='montgomery', algorithm='factored').codomain()
        for (c, isP) in [(c1b, isP1b), (c1f, isP1f), (c2b, isP2b), (c2f, isP2f)]:
            E = self._normalize_curve(E)
            P, Q = self._deterministic_torsion_basis(E, e)
            if isP:
                K = c * P + Q
            else:
                K = P + c * Q
            E = E.isogeny(K, model='montgomery', algorithm='factored').codomain()
        chl = self.Hash(msg + util.j_invariant_to_bytes(E) + util.j_invariant_to_bytes(Epk))
        return chl == c0 % 2**self.e_chl

    @staticmethod
    def _normalize_curve(E, points=()):
        """
        Algorithm 1 (MontgomeryNormalize) of the SQIsign specification
        (https://sqisign.org/spec/sqisign-20230601.pdf, Section 2.2.1.1).

        To a single supersingular j-invariant correspond six Montgomery
        A-invariants defining F_{p^2}-isomorphic curves. This returns the
        canonical one (depending only on the isomorphism class), together with
        the data (R, U) of an isomorphism

            E_A -> E_{A'},   (x, y) |-> (U^2 * (x + R), U^3 * y).
        """
        A = E.a2()
        F = A.parent()
        i = F.gen()
        A2 = A**2
        s = util.deterministic_sqrt(A2 - 4)
        u = (9 - A2) / 2
        t = (A**3 - 3*A) / (2 * s)
        # The three values Z0, Z1, Z2 are the squares of the six A-invariants.
        Z = min((A2, u + t, u - t), key=util.fp2_order_key)
        Aprime = util.deterministic_sqrt(Z)
        En = EllipticCurve(F, [0, Aprime, 0, 1, 0])
        if len(points) == 0:
            return En

        if Aprime == A:
            R, U = F(0), F(1)
        elif Aprime == -A:
            R, U = F(0), i
        else:
            R = (A2 + Aprime**2 - 6) * A / (A2 + 2*Aprime**2 - 9)
            U = util.deterministic_sqrt(Aprime / (A - 3*R))
        U2, U3 = U**2, U**3
        return En, [En(U2 * (T[0] + R), U3 * T[1]) for T in points]

    @staticmethod
    def _deterministic_torsion_basis(E, e):
        """
        Compute a deterministic basis of E[2^e] for a supersingular Montgomery curve.

        E is assumed to already be in canonical form (see MontgomeryNormalize): the
        caller must normalize the curve beforehand so that the basis depends only on
        the F_{p^2}-isomorphism class. No isomorphism is applied here.
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
