from sage.all import (
    sqrt,
    kronecker,
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
