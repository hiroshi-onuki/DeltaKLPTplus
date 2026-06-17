from sage.all import (
    sqrt,
    EllipticCurve,
)
import special_curve

class SQIsign:
    def __init__(self, p, e, f):
        E0withEnd = special_curve.SpecialSuperSingularCurve(p, e, f)
        self.E0withEnd = E0withEnd

    @staticmethod
    def _deterministic_torsion_basis(E, e):
        """
        Compute a deterministic basis of E[2^e] using Pari's implementation
        of the Weil pairing. This is used to compute the action of the
        quaternion order on E[2^e].
        """
        F = E.base_ring()
        p = F.characteristic()
        i = sqrt(F(-1))
        A = E.a2()
        assert A != 0
        assert E == EllipticCurve(F, [0, A, 0, 1, 0]) # y^2 = x^3 + A*x^2 + x

        h = F(0)
        xP = F(0)
        if A.is_square():
            print("A is square")
            while (1 + h**2).is_square() or not (xP**3 + A*xP**2 + xP).is_square():
                h += 1
                xP = -1/(1 + i*h)*A
        else:
            print("A is not square")
            while not (xP**3 + A*xP**2 + xP).is_square():
                h += 1
                xP = h*A
        cofactor = (p + 1) // (2**e)
        P = cofactor * E.lift_x(xP)
        Q = cofactor * E.lift_x(-xP - A)
        assert P.weil_pairing(Q, 2**e)**(2**(e - 1)) == -1
        return P, Q
