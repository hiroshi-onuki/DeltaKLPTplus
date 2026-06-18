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
)
from util import BiDLP_matrix_power_two
from utilities.discrete_log import tate_pairing_pari

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

    def KernelToIdeal(self, a, b):
        """
        return a left O-ideal I s.t. E0[I] = <a*P + b*Q>
        find a, b s.t.
            a*R + b*(qj + (1 + qk)/2)(R) = i(R), wehre R = a*P + b*Q
        """
        R = ZZ.quotient_ring(ZZ(2**self.e))
        M = 2*self.matrix_qi_qj - self.matrix_qi + self.matrix_1_qk # the action of (qj + (1 + qk)/2)
        M = matrix(R, M)
        v = vector(R, [a, b])
        Mv = M.transpose() * v
        M = matrix(R, [v, Mv]).transpose()
        v = M.inverse() * self.matrix_qi.transpose() * v
        a, b = [ZZ(c) for c in v]
        return self.order.left_ideal([a + b*(self.qj + (1 + self.qk)/2) - self.qi, 2**self.e])