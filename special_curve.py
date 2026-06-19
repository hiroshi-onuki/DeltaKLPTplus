from sage.all import (
    GF,
    ZZ,
    QuaternionAlgebra,
    EllipticCurve,
    PolynomialRing,
    is_square,
    matrix,
    identity_matrix,
    vector,
    inverse_mod,
    norm,
)
from util import BiDLP_matrix_power_two
from utilities.discrete_log import tate_pairing_pari
from theta_structures.couple_point import CouplePoint
from theta_isogenies.product_isogeny import EllipticProductIsogeny
from quaternion import Qlapoti

class SpecialSuperSingularCurve:
    def __init__(self, p, e, f):
        assert p == 2**e * f - 1
        assert e >= 3
        assert f % 2 == 1
        self.p = p
        self.e = e
        self.f = f
        Fpx = PolynomialRing(GF(p), 'x')
        x = Fpx.gen()
        self.Fp2 = GF(p**2, modulus=x**2 + 1, names='i')
        i = self.Fp2.gen()
        assert i**2 == -1
        E = EllipticCurve(self.Fp2, [1, 0]) # y^2 = x^3 + x
        self.curve = E
        B = QuaternionAlgebra(-1, -p)
        _, qi, qj, qk = B.basis()
        assert qi**2 == -1 and qj**2 == -p and qk == qi*qj
        O = B.quaternion_order([1, qi, (qi + qj)/2, (1 + qk)/2])
        assert O.is_maximal()
        self.order = O
        self.qi = qi
        self.qj = qj
        self.qk = qk

        """
        For computing the actions of qi, (qi + qj)/2, (1 + qk)/2 on E[2^e],
        we use a basis (Pext, Qext) of E[2^(e+1)] over GF(p^4)
        """
        tmp = self.Fp2.random_element()
        while is_square(tmp):
            tmp = self.Fp2.random_element()
        Fp2x = PolynomialRing(self.Fp2, 'x')
        x = Fp2x.gen()
        f = tmp.minpoly()(x**2)
        Fp4 = GF(p**4, modulus=f, names='j')
        Eext = EllipticCurve(Fp4, [1, 0])
        pi = Eext.frobenius_isogeny(1)
        Pext, Qext = Eext.torsion_basis(2**(e+1))

        # (P, Q) is a basis of E[2^e] over GF(p^2)
        emb = self.Fp2.embeddings(Fp4)[0]
        res = emb.section()
        def restrict_point(Pext):
            x, y = Pext.xy()
            return E([res(c) for c in (x, y)])
        P = restrict_point(2*Pext)
        Q = restrict_point(2*Qext)
        self.P = P
        self.Q = Q
        tPQ = tate_pairing_pari(P, Q, 2**e)**((p**2 - 1) // 2**e)
        self.tate_pairing_PQ = tPQ
        
        # The action of qi
        def qi_action(P):
            x, y = P.xy()
            return E([-x, i*y])
        iP = qi_action(P)
        iQ = qi_action(Q)
        self.matrix_qi = BiDLP_matrix_power_two(iP, iQ, P, Q, e)

        # The action of (qi + qj)/2
        def qi_action_ext(P):
            x, y = P.xy()
            return Eext([-x, emb(i)*y])
        Pd = restrict_point(qi_action_ext(Pext) + pi(Pext))
        Qd = restrict_point(qi_action_ext(Qext) + pi(Qext))
        self.matrix_qi_qj = BiDLP_matrix_power_two(Pd, Qd, P, Q, e)

        # The action of (1 + qk)/2
        Pd = restrict_point(Pext + qi_action_ext(pi(Pext)))
        Qd = restrict_point(Qext + qi_action_ext(pi(Qext)))
        self.matrix_1_qk = BiDLP_matrix_power_two(Pd, Qd, P, Q, e)

    def quaternion_action(self, alpha):
        a, b, c, d = vector(alpha) * self.order.basis_matrix().inverse()
        M = a * identity_matrix(2) + b * self.matrix_qi + c * self.matrix_qi_qj + d * self.matrix_1_qk
        a, b, c, d = M.list()
        return a * self.P + b * self.Q, c * self.P + d * self.Q

    def IdealToIsogeny(self, I):
        E0 = self.curve
        O = self.order
        N = norm(I)
        assert I.left_order() == O
        assert N % 2 == 1
        e = self.e

        """
        compute beta1, beta2, gamma s.t.
        beta1 = hat{phi_I} * phi_1, beta2 = hat{phi_I} * phi_2,
        gamma = hat{phi_2} * phi_1,
        where phi_I is the isogeny corresponding to I,
        deg(phi_1) + deg(phi_2) = 2^(e-2).
        """
        beta1, beta2, gamma = Qlapoti(I, e - 2)
        assert beta1 in I and beta2 in I
        assert beta2 * beta1.conjugate() / N == gamma
        d1 = beta1.reduced_norm() / N
        d2 = beta2.reduced_norm() / N
        assert d1 + d2 == 2**(e - 2)
        assert d1 % 2 == 1 and d2 % 2 == 1

        # compute (2^(e-2), 2^(e-2))-isogeny
        P1, Q1 = d1 * self.P, d1 * self.Q
        P2, Q2 = self.quaternion_action(gamma)
        K1 = CouplePoint(P1, P2)
        K2 = CouplePoint(Q1, Q2)
        Phi = EllipticProductIsogeny((K1, K2), e-2)
        P0, Q0 = self.quaternion_action(beta1.conjugate()) # P0 = hat{phi_1}*phi_I(P), Q0 = hat{phi_1}*phi_I(Q)
        image1 = Phi(CouplePoint(P0, E0(0)))
        image2 = Phi(CouplePoint(Q0, E0(0)))
        image_sum = Phi(CouplePoint(P0 + Q0, E0(0)))

        # determine the codomain EI by comparing the Tate pairings
        tP0Q0d1 = self.tate_pairing_PQ**(d1**2*N)
        for idx in range(2):
            Pim, Qim, PQim = image1[idx], image2[idx], image_sum[idx]
            if not (Pim + Qim == PQim or Pim + Qim == -PQim):
                Qim = -Qim
            assert Pim + Qim == PQim or Pim + Qim == -PQim
            exp = (self.p**2 - 1) / 2**e
            tPimQim = tate_pairing_pari(Pim, Qim, 2**e)**exp
            if tPimQim == tP0Q0d1:
                EI = Phi.codomain()[idx]
                PI = inverse_mod(d1, 2**e) * Pim
                QI = inverse_mod(d1, 2**e) * Qim
                return EI, PI, QI
        raise ValueError("Failed to determine the codomain of the isogeny")


    def KernelToIdeal(self, a, b, exp):
        """
        return a left O-ideal I s.t. E0[I] = <a*P' + b*Q'>, where (P', Q') = 2^(e-exp)*(P, Q)
        find a, b s.t.
            a*R + b*(qj + (1 + qk)/2)(R) = i(R), wehre R = a*P + b*Q
        """
        R = ZZ.quotient_ring(ZZ(2**exp))
        M = 2*self.matrix_qi_qj - self.matrix_qi + self.matrix_1_qk # the action of (qj + (1 + qk)/2)
        M = matrix(R, M)
        v = vector(R, [a, b])
        Mv = M.transpose() * v
        M = matrix(R, [v, Mv]).transpose()
        v = M.inverse() * self.matrix_qi.transpose() * v
        a, b = [ZZ(c) for c in v]
        return self.order.left_ideal([a + b*(self.qj + (1 + self.qk)/2) - self.qi, 2**exp])