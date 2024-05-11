def inversion(a):
    q = a.parent().order()
    r = bin(q-2)[3:]
    n = len(r)
    result = a
    for i in range(n):
        result = result ** 2
        if r[i] == '1':
            result *= a
    return result