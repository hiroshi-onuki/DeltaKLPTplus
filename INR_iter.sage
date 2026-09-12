
def INR_iter(p):
    N = p**5
    return float(4*log(p) * ceil(1/2 * log(log(RR(N/p) + 0.303, 2), 2) + 0.861))

