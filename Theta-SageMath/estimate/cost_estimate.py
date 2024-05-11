from time import process_time_ns
from sage.all import random_prime, randint, mod, GF, inverse_mod, set_random_seed

def cost_inverse(a):
    start = process_time_ns()
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
    bit_len = [254, 381, 1293] # the same with msi.rs
    print(f"{'Bit':>12}{'M':>12}{'S':>12}{'I':>12}")
    for bits in bit_len:
        lb = 1 << (bits-1)
        rb = (lb << 1) - 1
        p = random_prime(lbound=lb, n=rb)
        F = GF(p**2, name='i')
        i = F.gen()
        a = F(randint(lb, rb) * i + randint(lb, rb))
        b = F(randint(lb, rb) * i + randint(lb, rb))

        num = 100
        sum_multiple = 0
        sum_square = 0
        sum_inverse = 0
        for _ in range(num):
            sum_multiple += cost_multiple(a, b)
            sum_square += cost_square(a)
            sum_inverse += cost_inverse(a)
            # sum_inverse += cost_inverse_old(a)
        print(f"{bits:>12}{sum_multiple/num:>12}{sum_square/num:>12}{sum_inverse/num:>12}")
