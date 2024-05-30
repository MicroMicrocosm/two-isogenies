from time import process_time_ns
from sage.all import randint, GF, ZZ, set_random_seed

from isogeny_diamond import DIAMONDS

def cost_inverse(a):
    start = process_time_ns()
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

    end = process_time_ns()
    return end - start

def cost_inverse_old(a):
    start = process_time_ns()
    n = 1 / a
    end = process_time_ns()
    return end - start

def cost_square(a):
    start = process_time_ns()
    n = a * a
    end = process_time_ns()
    return end - start

def cost_multiple(a, b):
    start = process_time_ns()
    n = a * b
    end = process_time_ns()
    return end - start

if __name__ == "__main__":
    set_random_seed(0)
    # bit_len = [254, 381, 1293] # the same with msi.rs
    print(f"{'Bit':>12}{'M':>12}{'S':>12}{'S/M':>12}{'I':>12}{'I/M':>12}{'I(old)':>12}{'I/M(old)':>12}")
    for f, ea, eb, _, _ in DIAMONDS:
        A = ZZ(2 ** ea)
        B = ZZ(3 ** eb)
        p = 4 * f * A * B - 1
        bit = p.nbits()
        F = GF(p**2, name='i', modulus=[1, 0, 1])
        i = F.gen()
        a = F(randint(0, p-1) * i + randint(0, p-1))
        b = F(randint(0, p-1) * i + randint(0, p-1))

        num = 100
        sum_multiple = 0
        sum_square = 0
        sum_inverse = 0
        sum_inverse_old = 0
        for _ in range(num):
            sum_multiple += cost_multiple(a, b)
            sum_square += cost_square(a)
            sum_inverse += cost_inverse(a)
            sum_inverse_old += cost_inverse_old(a)
        print(f"{bit:>12}{sum_multiple/num:>12}{sum_square/num:>12}{sum_square/sum_multiple:>12.3f}{sum_inverse/num:>12}{sum_inverse/sum_multiple:>12.3f}{sum_inverse_old/num:12}{sum_inverse_old/sum_multiple:>12.3f}")
