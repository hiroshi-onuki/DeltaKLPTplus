from sage.all import (
    matrix,
    ZZ,
)
from utilities.discrete_log import BiDLP_power_two

def BiDLP_matrix_power_two(R, S, P, Q, e):
    """
    return the matrix [[a, b], [c, d]] such that
    R = a*P + b*Q, S = c*P + d*Q, where P, Q is a basis of E[2^e].
    """
    a, b = BiDLP_power_two(R, P, Q, e, None)
    c, d = BiDLP_power_two(S, P, Q, e, None)
    return matrix(ZZ, 2, 2, [a, b, c, d])
