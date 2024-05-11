from sage.all import ZZ
from utilities.discrete_log import weil_pairing_pari
from utilities.inversion import inversion
from montgomery_isogenies.kummer_line import KummerLine, KummerPoint


class CouplePoint:
    """
    A helper class which represents an element P = (P1, P2) in E1 x E2
    and allows us to compute certain useful functions, such as adding,
    doubling or comouting the Weil pairing of e(P,Q) for P,Q in E1 x E2
    """

    def __init__(self, P1, P2):
        self.P1 = P1
        self.P2 = P2

    def __repr__(self):
        return "[{},{}]".format(self.P1, self.P2)

    def parent(self):
        return (self.P1.curve(), self.P2.curve())

    def curves(self):
        return self.parent()

    def points(self):
        return self.P1, self.P2

    def order(self):
        return (self.P1.order(), self.P2.order())

    def double(self):
        """
        Computes [2] P = ([2] P1, [2] P2)
        """
        # TODO
        return ZZ(2) * self

    def double_iter_old(self, n):
        """
        Compute [2^n] P = ([2^n] P1, [2^n] P2)
        """
        # When the scalar is a python int, then
        # sagemath does multiplication naively, when
        # the scalar in a Sage type, it instead calls
        # _acted_upon_, which calls pari, which is fast
        m = ZZ(2**n)
        return m * self
    
    def double_iter(self, n):
        """
        Compute [2^n] P = ([2^n] P1, [2^n] P2)
        
        Using x-only coordinates and Okeya-Sakurai formula when computing [2^n] Pi
        """
        P_2n = ()
        for P in self.points():
            curve = P.curve()
            base_ring = curve.base_ring()

            ainvs = curve.a_invariants()
            A = ainvs[1]
            if ainvs != (0, A, 0, 1, 0):
                raise ValueError("Must be Montgomery curve.")
            A = base_ring(A)
            
            P_Kummer = KummerPoint(KummerLine(curve), P)
            X1, Z1, X0, Z0 = P_Kummer.double_iter(n) # cost : 6M + 4S + 1C
            xP, yP = P[0], P[1]

            # recover (x1, y1) cost : 12M + 3S + 1I
            y1 = X1 * xP + Z1
            t0 = X1 + xP * Z1
            y1 = t0 * y1
            y1 = y1 + y1
            t0 = t0 ** 2
            t1 = X1 - xP * Z1
            t1 = t1 ** 2
            t0 = t0 - t1
            t0 = A * t0
            y1 = y1 + t0
            y1 = Z0 * y1
            t1 = X0 * t1
            t1 = t1 + t1
            y1 = t1 - y1
            t0 = yP * Z0
            x1 = X1 * Z1
            x1 = t0 * x1
            x1 = x1 + x1
            x1 = x1 + x1
            t1 = Z1 ** 2
            t0 = t0 * t1
            t0 = t0 + t0
            t0 = t0 + t0
            t0 = inversion(t0)
            x1 = x1 * t0
            y1 = y1 * t0

            P_2n = P_2n + (curve(x1, y1), )
        
        return CouplePoint(P_2n[0], P_2n[1])

    def __getitem__(self, i):
        # Operator to get self[i].
        if i == 0:
            return self.P1
        elif i == 1:
            return self.P2
        else:
            raise ValueError("Index {} is out of range.".format(i))

    def __setitem__(self, i, P):
        # Operator to set self[i]=P.
        if i == 0:
            self.P1 = P
        elif i == 1:
            self.P2 = P
        else:
            raise ValueError("Index {} is out of range.".format(i))

    def __eq__(self, other):
        return self.P1 == other.P1 and self.P2 == other.P2

    def __add__(self, other):
        return CouplePoint(self.P1 + other.P1, self.P2 + other.P2)

    def __sub__(self, other):
        return CouplePoint(self.P1 - other.P1, self.P2 - other.P2)

    def __neg__(self):
        return CouplePoint(-self.P1, -self.P2)

    def __mul__(self, m):
        """
        Compute [m] P = ([m] P1, [m] P2)
        """
        # When the scalar is a python int, then
        # sagemath does multiplication naively, when
        # the scalar in a Sage type, it instead calls
        # _acted_upon_, which calls pari, which is fast
        m = ZZ(m)
        return CouplePoint(m * self.P1, m * self.P2)

    def __rmul__(self, m):
        return self * m

    def weil_pairing(self, other, n):
        """
        The Weil pairing e_n(P, Q) for P = (P1, P2) and Q = (Q1, Q2)
        is defined as

            e_n(P, Q) = e_n(P1, Q1) * e_n(P2, Q2)
        """
        if not isinstance(other, CouplePoint):
            raise TypeError("Both inputs must be couple points")

        P1, P2 = self.points()
        Q1, Q2 = other.points()

        ePQ1 = weil_pairing_pari(P1, Q1, n)
        ePQ2 = weil_pairing_pari(P2, Q2, n)

        Fp2 = P1.base_ring()
        return Fp2(ePQ1 * ePQ2)
