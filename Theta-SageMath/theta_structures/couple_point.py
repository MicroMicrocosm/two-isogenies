from sage.all import ZZ
from utilities.discrete_log import weil_pairing_pari
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
    
    @staticmethod
    def x_dbl_add(XP, ZP, XQ, ZQ, xPQ, zPQ, a):
        """
        function for step in Montgomery ladder
        simultaneous doubling and differential addition
        
        Input: projective coordinates P=(XP:ZP) and Q=(XQ:ZQ), 
               projective difference P-Q=(xPQ:zPQ) and 
               curve constant a = (A+2)/4.   
        Output: projective coordinates of 2P=(X2P:Z2P)
                and Q+P=(XQP:ZQP)

        Cost: 4S + 6M + 1C
        """
        
        t0 = XP + ZP                  
        t1 = XP - ZP 
        X2P = t0 * t0
        t2 = XQ - ZQ
        XQP = XQ + ZQ
        t0 = t0 * t2
        Z2P = t1 * t1
        t1 = t1 * XQP
        t2 = X2P - Z2P
        X2P = X2P * Z2P
        XQP = a * t2
        ZQP = t0 - t1
        Z2P = XQP + Z2P
        XQP = t0 + t1
        Z2P = Z2P * t2
        ZQP = ZQP * ZQP
        XQP = XQP * XQP
        ZQP = xPQ * ZQP
        XQP = XQP * zPQ

        return X2P, Z2P, XQP, ZQP

    def double_iter(self, n, flag=False):
        """
        Compute [2^n] P = ([2^n] P1, [2^n] P2)
        
        Using x-only coordinates and Okeya-Sakurai formula when computing [2^n] Pi deponds on whether flag is True.

        Cost : 8S + 12M + 2C if flag is True else 12S + 10M + 2C
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

            if flag:
                a = (A + 2) / 4
                XP, ZP = P[0], P[2]
                X2n, Z2n = P[0], P[2]
                X2n1, Z2n1 = base_ring.one(), base_ring.zero()
                for _ in range(n):
                    X2n, Z2n, X2n1, Z2n1 = self.x_dbl_add(X2n, Z2n, X2n1, Z2n1, XP, ZP, a)
                x1, y1 = P[0], P[1]

                # recover [2^n]P = (X : Y : Z) cost : 3S + 11M
                t = Z2n * Z2n1
                X = X2n * t
                Z = Z2n * t
                t = x1 * x1 + 1
                Y = A * x1
                Y = Y + Y
                t = t + Y
                t = t * X
                Y = X2n * X2n
                Y = Y * Z2n1
                Y = Y + Z
                Y = Y * x1
                Y = Y + t
                t = x1 * Z2n
                t = X2n - t
                t = t * t
                t = t * X2n1
                Y = t - Y
                t = y1 + y1
                X = X * t
                Z = Z * t
            else:
                X, Y, Z = P
                for _ in range(n):
                    xx = X * X
                    zz = Z * Z
                    dxz = (X + Z) * (X + Z) - xx - zz
                    dyz = 2 * Y * Z
                    t0 = xx - zz
                    t1 = xx + zz
                    X = t0 * t0
                    X = dyz * X
                    Y = t1 + A * dxz
                    Y = t1 * Y
                    Y = Y + dxz * dxz
                    Y = t0 * Y
                    Z = dyz * dyz
                    Z = dyz * Z

            P_2n = P_2n + (curve(X, Y, Z), )
        
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
