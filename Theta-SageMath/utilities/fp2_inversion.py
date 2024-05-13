def inversion(a):
    p = a.parent().characteristic()
    r = bin(p-2)[3:]
    n = len(r)

    i = a.parent().gen()
    a_conj = a.conjugate()  # a_conj = x0 - y0*i
    x = a + a_conj          # x = 2 * x0
    y = i * (a_conj - a)    # y = 2 * y0
    
    b = x * x + y * y
    result = b
    for k in range(n):
        result = result * result
        if r[k] == '1':
            result *= b
    result = result * a_conj
    result = result + result
    result = result + result
    return result