def invert_Fp2(a):
    p = a.parent().characteristic()
    r = bin(p-2)[3:]
    n = len(r)

    x, y = a.list()             # a = x + yi
    if y: # y != 0
        a_conj = a.conjugate()  # a_conj = x - yi
        b = x * x + y * y
    else:
        a_conj = 1
        b = x

    result = b
    for k in range(n):
        result = result * result
        if r[k] == '1':
            result *= b
    result = result * a_conj
    return result