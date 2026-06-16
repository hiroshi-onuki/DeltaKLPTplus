from sage.all import (
    EllipticCurve,
    GF,
    PolynomialRing,
    is_square,
    discrete_log,
)
from utilities import discrete_log

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
        self.E = EllipticCurve(self.Fp2, [1, 0]) # y^2 = x^3 + x

        # compute the actions of qi, (qi + qj)/2, (1 + qk)/2 on E[2^e]
        tmp = self.Fp2.random_element()
        while is_square(tmp):
            tmp = self.Fp2.random_element()
        Fp2x = PolynomialRing(self.Fp2, 'x')
        x = Fp2x.gen()
        Fp4 = self.Fp2.extension(x**2 - tmp)
        Ext = self.E.base_extend(Fp4)
        P, Q = Ext.torsion_basis(2**(e+1))
        print(f"Found basis P, Q of E[2^{e+1}] over Fp4")

        R = 10*P + 15*Q
        a, b = self._BiDLP(R, P, Q, 2**(e+1))
        print(f"Action of qi: {a}P + {b}Q")

    def _BiDLP(R, P, Q, n):
        # solve Bi-DLP R = aP + bQ
        ePQ = P.weil_pairing(Q, n)
        ePR = P.weil_pairing(R, n)
        eRQ = R.weil_pairing(Q, n)

        print(f"e(P, Q) = {ePQ}, e(P, R) = {ePR}, e(R, Q) = {eRQ}")
        a = discrete_log(eRQ, ePQ)
        b = discrete_log(ePR, ePQ)
        print(f"Found a={a}, b={b} such that e(R, Q) = e(P, Q)^a and e(P, R) = e(P, Q)^b")
        assert R == a*P + b*Q
        return a, b

