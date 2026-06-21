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

def integer_to_bytes(n, byte_len=None):
    """
    Represent an integer n as bytes in big-endian
    """
    n = int(n)
    if byte_len is None:
        byte_len = (n.bit_length() + 7) // 8
    return n.to_bytes(byte_len, 'big')

def bytes_to_integer(b):
    """
    Represent bytes as an integer in big-endian
    """
    return ZZ(int.from_bytes(b, 'big'))

def field_element_to_bytes(x, byte_len):
    """
    Represent an element x = a + i*b in Fp^2 as bytes
    """
    # x = a + b*i
    a_int, b_int = x.list()
    # convert to bytes
    a_bytes = integer_to_bytes(a_int, byte_len=byte_len)
    b_bytes = integer_to_bytes(b_int, byte_len=byte_len)
    # concatenate bytes
    return a_bytes + b_bytes

def field_element_from_bytes(Fp2, x_bytes, byte_len):
    """
    Recover an element x = a + i*b in Fp^2 from bytes
    """
    a_bytes = x_bytes[:byte_len]
    b_bytes = x_bytes[byte_len:]
    a_int = bytes_to_integer(a_bytes)
    b_int = bytes_to_integer(b_bytes)
    return Fp2([a_int, b_int])

def j_invariant_to_bytes(E):
    """
    Represent the j-invariant of an elliptic curve E as bytes
    """
    p = E.base_ring().characteristic()
    byte_len = (len(bin(p)) - 2 + 7) // 8
    j = E.j_invariant()
    return field_element_to_bytes(j, byte_len)