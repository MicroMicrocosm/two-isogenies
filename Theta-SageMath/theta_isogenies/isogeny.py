from sage.all import ZZ

from theta_structures.dimension_two import ThetaStructure, ThetaPoint
from theta_isogenies.morphism import Morphism
from utilities.batched_inversion import batched_inversion


class ThetaIsogeny(Morphism):
    def __init__(self, domain, T1_8, T2_8, hadamard=(False, True), flag=False):
        """
        Compute a (2,2)-isogeny in the theta model. Expects as input:

        - domain: the ThetaStructure from which we compute the isogeny
        - (T1_8, T2_8): points of 8-torsion above the kernel generating the isogeny

        When the 8-torsion is not available (for example at the end of a long
        (2,2)-isogeny chain), the the helper functions in isogeny_sqrt.py
        must be used.

        NOTE: on the hadamard bools:

        The optional parameter 'hadamard' controls if we are in standard or dual
        coordinates, and if the codomain is in standard or dual coordinates. By
        default this is (False, True), meaning we use standard coordinates on
        the domain A and the codomain B.

        The kernel is then the kernel K_2 where the action is by sign. Other
        possibilities: - (False, False): standard coordinates on A, dual
        coordinates on B - (True, True): start in dual coordinates on A
        (alternatively: standard coordinates on A but quotient by K_1 whose
        action is by permutation), and standard coordinates on B. - (True,
        False): dual coordinates on A and B

        These can be composed as follows for A -> B -> C:

        - (False, True) -> (False, True) (False, False) -> (True, True):
          - standard coordinates on A and C,
          - standard/resp dual coordinates on B
        - (False, True) -> (False, False) (False, False) -> (True, False):
          - standard coordinates on A,
          - dual coordinates on C,
          - standard/resp dual coordinates on B
        - (True, True) -> (False, True) (True, False) -> (True, True):
          - dual coordinates on A,
          - standard coordinates on C,
          - standard/resp dual coordiantes on B
        - (True, True) -> (False, False) (True, False) -> (True, False):
          - dual coordinates on A and C
          - standard/resp dual coordinates on B

        On the other hand, these gives the multiplication by [2] on A:

        - (False, False) -> (False, True) (False, True) -> (True, True):
          - doubling in standard coordinates on A
          - going through dual/standard coordinates on B=A/K_2
        - (True, False) -> (False, False) (True, True) -> (True, False):
          - doubling in dual coordinates on A
          - going through dual/standard coordinates on B=A/K_2
            (alternatively: doubling in standard coordinates on A going
            through B'=A/K_1)
        - (False, False) -> (False, False) (False, True) -> (True, False):
          - doubling from standard to dual coordinates on A
        - (True, False) -> (False, True) (True, True) -> (True, True):
          - doubling from dual to standard coordinates on A
        """
        if not isinstance(domain, ThetaStructure):
            raise ValueError
        self._domain = domain

        self._hadamard = hadamard
        self.flag = flag
        self._precomputation = None
        self._codomain = self._compute_codomain(T1_8, T2_8)

    def _compute_codomain(self, T1, T2):
        """
        Given two isotropic points 8-torsion T1 and T2, compatible with
        the theta null point, compute the level two theta null point A/K_2

        Cost : 8S + 9M if flag is True else 8S + 13M + 1I / 8S + 23M + 1I without precomputation
        """
        if self._hadamard[0]:
            xA, xB, _, _ = ThetaPoint.to_squared_theta(
                *ThetaPoint.to_hadamard(*T1.coords())
            )
            zA, tB, zC, tD = ThetaPoint.to_squared_theta(
                *ThetaPoint.to_hadamard(*T2.coords())
            )
        else:
            xA, xB, _, _ = T1.squared_theta()
            zA, tB, zC, tD = T2.squared_theta()

        if self.flag:
            # Compute A, B, C, D
            xAtB = xA * tB
            zAxB = zA * xB
            A = xAtB * zA
            B = zAxB * tB
            C = xAtB * zC
            D = zAxB * tD

            # Compute the inverse of A, B, C, D
            zCtD = zC * tD
            A_inv = xB * zCtD
            B_inv = xA * zCtD
            C_inv = D
            D_inv = C
        else:
            if not self._hadamard[0] and self._domain._precomputation:
                # Batch invert denominators
                xA_inv, zA_inv, tB_inv = batched_inversion(xA, zA, tB)

                # Compute A, B, C, D
                A = ZZ(1)
                B = xB * xA_inv
                C = zC * zA_inv
                D = tD * tB_inv * B

                _, _, _, _, _, BBinv, CCinv, DDinv = self._domain.precomputation(self.flag)
                A_inv = ZZ(1)
                B_inv = BBinv * B
                C_inv = CCinv * C
                D_inv = DDinv * D
            else:
                # Batch invert denominators
                xA_inv, zA_inv, tB_inv, xB_inv, zC_inv, tD_inv = batched_inversion(
                    xA, zA, tB, xB, zC, tD
                )

                # Compute A, B, C, D
                A = ZZ(1)
                B = xB * xA_inv
                C = zC * zA_inv
                D = tD * tB_inv * B

                A_inv = ZZ(1)
                B_inv = xB_inv * xA
                C_inv = zC_inv * zA
                D_inv = tD_inv * tB * B_inv

        self._precomputation = (A_inv, B_inv, C_inv, D_inv)
        if self._hadamard[1]:
            a, b, c, d = ThetaPoint.to_hadamard(A, B, C, D)
            return ThetaStructure([a, b, c, d])
        else:
            return ThetaStructure([A, B, C, D])
    
    def __call__(self, P):
        """
        Take into inout the theta null point of A/K_2, and return the image
        of the point by the isogeny

        Cost : 4S + 4M if flag is True else 4S + 3M
        """
        if not isinstance(P, ThetaPoint):
            raise TypeError("Isogeny evaluation expects a ThetaPoint as input")
        
        if self._hadamard[0]:
            xx, yy, zz, tt = ThetaPoint.to_squared_theta(
                *ThetaPoint.to_hadamard(*P.coords())
            )
        else:
            xx, yy, zz, tt = P.squared_theta()

        if self.flag:
            A_inv, B_inv, C_inv, D_inv = self._precomputation
            xx = xx * A_inv
            yy = yy * B_inv
            zz = zz * C_inv
            tt = tt * D_inv
        else:
            _, B_inv, C_inv, D_inv = self._precomputation
            yy = yy * B_inv
            zz = zz * C_inv
            tt = tt * D_inv

        image_coods = (xx, yy, zz, tt)
        if self._hadamard[1]:
            image_coods = ThetaPoint.to_hadamard(*image_coods)
        return self._codomain(image_coods)
