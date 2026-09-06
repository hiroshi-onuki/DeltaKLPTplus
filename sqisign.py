from sage.all import ( # type: ignore
    ZZ,
    is_prime,
    sqrt,
    kronecker,
    EllipticCurve,
    vector,
    norm,
    inverse_mod,
    randint,
    log,
    ceil,
    pi,
)
import hashlib
import special_curve
import quaternion
import util
import montgomery
import utilities.discrete_log

class SQIsign:
    def __init__(self, e, f, lam):
        self.sec_lambda = ZZ(lam)
        self.e = ZZ(e)
        self.f = ZZ(f)
        p = 2**self.e * self.f - 1
        assert is_prime(p)
        self.p = p
        self.E0withEnd = special_curve.SpecialSuperSingularCurve(p, self.e, self.f)
        Dmix = p * 2**(2*self.sec_lambda) + 1
        while not is_prime(Dmix):
            Dmix += 2
        self.Dmix = Dmix    # the degree of phi_sk and phi_com, which satisfies the mixing property in the supersingular isogeny graph
        self.e_chl = self.sec_lambda
        self.e_rsp = self._response_length(p, self.sec_lambda)

    def Keygen(self):
        Isk, _ = quaternion.RandomFixedNormIdeal(self.E0withEnd.order, self.Dmix)
        Epk, Psk, Qsk = self.E0withEnd.IdealToIsogeny(Isk)
        Epk, (Psk, Qsk) = self._normalize_curve(Epk, (Psk, Qsk))
        Ppk, Qpk = self._deterministic_torsion_basis(Epk, self.E0withEnd.e)
        Msk = util.BiDLP_matrix_power_two(Psk, Qsk, Ppk, Qpk, self.E0withEnd.e)
        sk = (Isk, Msk)
        pk = Epk
        return sk, pk

    def Commit(self):
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

        IcomIrsp, _ = quaternion.DeltaKLPT_plus(Icom, IskIchl, 2, self.e_rsp, self.sec_lambda)
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

        # the coefficient for the challenge + the remaining part of the response
        e0 = self.e_chl + self.e_rsp - 4*e
        Im0p2 = Iall + O0 * 2**e0
        v0 = self.E0withEnd.IdealToKernel(Im0p2, e0)
        v0 = v0 * Msk % 2**e0
        c0 = v0[1] * inverse_mod(v0[0], 2**e0) % 2**e0
        assert c0 % 2**self.e_chl == chl
        c0without_chl = (c0 - chl) // 2**self.e_chl

        # the coefficients for 2^(4e)-part of the response
        nIsk = norm(Isk)
        coords = []
        for k in (1, 2):
            Im, Imp2b, Imp2f = self._setup_response_ideal(Iall, k, e0, e, nIsk, O0)
            coords.append(self._response_coords(Im, Imp2b, Imp2f, e))
        (c1b, c1f, isP1b, isP1f), (c2b, c2f, isP2b, isP2f) = coords

        return (c0without_chl, c1b, c1f, c2b, c2f, isP1b, isP1f, isP2b, isP2f)

    def _setup_response_ideal(self, Iall, k, e0, e, nIsk, O0):
        """
        Build the k-th (k = 1, 2) response ideal Im together with the forward and
        backward 2^e-torsion push-out ideals (Imp2f, Imp2b) used to read off its
        coordinates. The two responses differ only by the exponent offset 2*k*e.
        """
        hi = 2**(e0 + 2*k*e) * nIsk
        lo = 2**(e0 + (2*k - 1)*e) * nIsk
        ImImp2f = Iall + O0 * hi
        Im = ImImp2f + O0 * lo
        Im, beta, _ = quaternion.EquivalentPrimeIdeal(Im, self.sec_lambda)
        Imp2f = ImImp2f * (beta.conjugate() / lo) + O0 * 2**e
        Imp2b = O0 * beta.conjugate() + O0 * 2**e
        return Im, Imp2b, Imp2f

    def _response_coords(self, Im, Imp2b, Imp2f, e):
        """
        Translate one response ideal into its (backward, forward) coordinates.

        Returns (cb, cf, isPb, isPf): cb/isPb describe the kernel pushed through
        the dual isogeny, cf/isPf the forward 2^e-torsion kernel.
        """
        Em, Pm, Qm = self.E0withEnd.IdealToIsogeny(Im)
        Em, (Pm, Qm) = self._normalize_curve(Em, (Pm, Qm))
        Pmd, Qmd = self._deterministic_torsion_basis(Em, e)
        Mm = util.BiDLP_matrix_power_two(Pm, Qm, Pmd, Qmd, e)

        # backward: evaluate the complementary point through the dual isogeny
        vdual = self.E0withEnd.IdealToKernel(Imp2b, e) * Mm % 2**e
        Kdual = vdual[0] * Pmd + vdual[1] * Qmd
        evalP = Pmd if vdual[0] % 2 == 0 else Qmd
        phi = Em.isogeny(Kdual, model='montgomery', algorithm='factored')
        K = phi(evalP)
        assert K.order() == 2**e
        Emb, (K,) = self._normalize_curve(phi.codomain(), (K,))
        Pmb, Qmb = self._deterministic_torsion_basis(Emb, e)
        a, b = utilities.discrete_log.BiDLP_power_two(K, Pmb, Qmb, e, None)
        cb, isPb = self._projective_coord(a, b, e)

        # forward: read the 2^e-torsion kernel directly
        vf = self.E0withEnd.IdealToKernel(Imp2f, e) * Mm % 2**e
        cf, isPf = self._projective_coord(vf[0], vf[1], e)

        return cb, cf, isPb, isPf

    @staticmethod
    def _projective_coord(a, b, e):
        """
        Encode the kernel spanned by (a, b) with a single coordinate.

        Returns (c, isP): with isP True the kernel is c*P + Q, otherwise P + c*Q,
        where c lives in Z/2^e. The even one of a, b is placed in the numerator so
        the other is invertible mod 2^e.
        """
        if a % 2 == 0:
            return a * inverse_mod(b, 2**e) % 2**e, True
        return b * inverse_mod(a, 2**e) % 2**e, False

    def RecoverCommitment(self, pk, chl, rsp):
        """
        Recompute the commitment curve from (pk, chl, rsp) by walking the five 2-power isogenies
        of the response, and report whether the whole walk is cyclic.

        The walk is done on the Kummer line (montgomery.py): each kernel <P + c*Q> or <c*P + Q>
        is obtained with a 3-point ladder from x(P), x(Q), x(P - Q) of the deterministic basis,
        the chain of 2-isogenies returns the codomain coefficient A and the x-only image of the
        cyclicity witness 2^(e-1) Q, and montgomery.normalize_A plays the role of
        _normalize_curve. Sage curves are only built where the deterministic basis is computed.
        """
        e = self.E0withEnd.e
        e0 = self.e_chl + self.e_rsp - 4*e
        c0without_chl, c1b, c1f, c2b, c2f, isP1b, isP1f, isP2b, isP2f = rsp
        c0 = c0without_chl * 2**self.e_chl + chl
        Epk = pk
        Epk.set_order((self.p + 1)**2, check=False)  # pk may come from outside this session
        A = Epk.a2()
        xP, xQ, xPQ = self._basis_x(Epk, e)
        A24p, C24 = montgomery.A24_projective(A)
        xK = montgomery.ladder3pt(c0, xP, xQ, xPQ, montgomery.A24_affine(A), nbits=e0)  # x(P + c0 Q)
        xK = montgomery.xDBLe(xK, A24p, C24, e - e0)                                      # order 2^e0
        imP = montgomery.xDBLe((xQ, 1), A24p, C24, e - 1)   # 2^(e-1) Q, for the cyclicity check
        A, (imP,) = montgomery.isogeny_chain_2e(A, xK, e0, [imP])
        for (c, isP) in [(c1b, isP1b), (c1f, isP1f), (c2b, isP2b), (c2f, isP2f)]:
            Aprime, R, U2 = montgomery.normalize_A(A)
            imP = montgomery.apply_isomorphism_x(imP, R, U2)
            En = self._curve_from_A(Aprime)
            xP, xQ, xPQ = self._basis_x(En, e)
            xK = self._kernel_x(c, isP, xP, xQ, xPQ, montgomery.A24_affine(Aprime), e)
            A, (imP,) = montgomery.isogeny_chain_2e(Aprime, xK, e, [imP])
        E = self._curve_from_A(A)
        return E, imP[1] != 0

    def _curve_from_A(self, A):
        """The Sage curve y^2 = x^3 + A x^2 + x with its order set."""
        F = self.E0withEnd.Fp2
        E = EllipticCurve(F, [0, A, 0, 1, 0])
        E.set_order((self.p + 1)**2, check=False)
        return E

    def _basis_x(self, E, e):
        """x(P), x(Q), x(P - Q) of the deterministic basis of E[2^e]."""
        P, Q = self._deterministic_torsion_basis(E, e)
        return P[0], Q[0], (P - Q)[0]

    def _kernel_x(self, c, isP, xP, xQ, xPQ, A24, e):
        """x(c*P + Q) if isP else x(P + c*Q), via the 3-point ladder (c < 2^e)."""
        if isP:
            return montgomery.ladder3pt(c, xQ, xP, xPQ, A24, nbits=e)
        return montgomery.ladder3pt(c, xP, xQ, xPQ, A24, nbits=e)

    def Hash(self, msg):
        h = hashlib.sha256(msg)
        c = util.bytes_to_integer(h.digest())
        return c % 2**self.e_chl

    def Sign(self, sk, pk, msg):
        com, st = self.Commit()
        chl = self.Hash(msg + util.j_invariant_to_bytes(com) + util.j_invariant_to_bytes(pk))
        rsp = self.Respond(sk, st, chl)
        c0without_chl = rsp[0]
        c0 = c0without_chl * 2**self.e_chl + chl
        rsp = (c0,) + rsp[1:]
        return rsp

    def Verify(self, pk, msg, sign):
        c0 = sign[0]
        chl = c0 % 2**self.e_chl
        # RecoverCommitment expects c0 without the challenge part (it re-adds chl);
        # the signature stores the full c0, so strip the low e_chl bits here.
        c0without_chl = c0 // 2**self.e_chl
        rsp = (c0without_chl,) + sign[1:]
        com, is_cyclic = self.RecoverCommitment(pk, chl, rsp)
        chl_check = self.Hash(msg + util.j_invariant_to_bytes(com) + util.j_invariant_to_bytes(pk))
        return chl == chl_check and is_cyclic
         
    def Simulator(self, pk, chl):
        e = self.E0withEnd.e
        e0 = self.e_chl + self.e_rsp - 4*e
        Epk = pk
        Epk.set_order((self.p + 1)**2, check=False)  # pk may come from outside this session
        Ppk, Qpk = self._deterministic_torsion_basis(Epk, self.E0withEnd.e)

        c0without_chl = randint(0, 2**(e0 - self.e_chl) - 1)
        c0 = c0without_chl * 2**self.e_chl + chl
        K = 2**(e - e0) * (Ppk + c0 * Qpk)
        phi = Epk.isogeny(K, model='montgomery', algorithm='factored')
        E = phi.codomain()
        imP = phi(2**(e - 1) * Qpk)  # for checking the cyclicity of the isogeny
        
        rsp_c = (c0without_chl,)
        rsp_isP = ()
        for _ in range(4):
            E, (imP,) = self._normalize_curve(E, (imP,))
            P, Q = self._deterministic_torsion_basis(E, e)
            P2, Q2 = 2**(e - 1) * P, 2**(e - 1) * Q
            c = randint(0, 2**e - 1)
            if imP == P2:
                if c % 2 == 0:
                    K = c * P + Q
                    isP = True
                else:
                    K = P + c * Q
                    isP = False
            elif imP == Q2:
                K = P + c * Q
                isP = False
            else:   # imP == P2 + Q2
                if c % 2 == 0:
                    K = c * P + Q
                    isP = True
                else:
                    c -= 1
                    K = P + c * Q
                    isP = False
            rsp_c += (c,)
            rsp_isP += (isP,)
            phi = E.isogeny(K, model='montgomery', algorithm='factored')
            E = phi.codomain()
            imP = phi(imP)
        assert not imP.is_zero()
        rsp = rsp_c + rsp_isP
        return E, rsp

    @staticmethod
    def _response_length(p, omega):
        return ceil(log(25*log(2)/(6*pi) * omega * p**4 * log(p), 2))

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
        # All curves in the scheme are supersingular with (p+1)^2 rational
        # points. Declaring the order here (isogeny codomains inherit it) avoids
        # an expensive point counting inside every point.order() call.
        En.set_order((F.characteristic() + 1)**2, check=False)
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
        return En, [En(0) if T.is_zero() else En(U2 * (T[0] + R), U3 * T[1]) for T in points]

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
