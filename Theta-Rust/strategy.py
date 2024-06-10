# ==================================================== #
#       Compute optimised strategy for (2,2)-chain     #
# ---------------------------------------------------- #
# Here we deviate from the usual formula by allowing   #
# doubling and images from the elliptic product to     #
# have a distinct cost from the rest of the tree,      #
# helping limit the number of costly images from the   #
# product in favour for slightly more doubling         #
# ==================================================== #

import functools
import sys

# For the long chain, we need a lot of recursion!
sys.setrecursionlimit(1500)

# fmt: off
def optimised_strategy_previous(n, M, S, I):
    """
    A modification of

    Algorithm 60: https://sike.org/files/SIDH-spec.pdf Shown to be appropriate
    for (l,l)-chains in https://ia.cr/2023/508

    Which allows the leftmost branch to have a different cost for the rest of
    the tree. This is partiularly useful for (2,2) isogenies, where the gluing
    doubling and images have a much higher cost than the rest of the tree.

    Thanks to Robin Jadoul for helping with the implementation of this function 
    via personal communication
    """

    # Define the costs and initalise the nodes which we store during doubling
    left_cost = (8*S + 6*M, 12*M + 12*S)       # (regular_cost, left_branch_cost) Double
    right_cost = (4*S + 3*M, 82*M + 18*S + I)  # (regular_cost, first_right_cost) Images
    checkpoints = ({}, {})  # (inner, left edge)

    @functools.cache
    def cost(n, leftmost):
        """
        The minimal cost to get to all children of a height `n` tree.
        If `leftmost` is true, we're still on the leftmost edge of the "outermost" tree

        Updates a global "Check points" which are the points along a branch which we 
        keep for later
        """
        if n <= 1:
            return 0  # no cost here

        c = float("inf")
        for i in range(1, n):  # where to branch off
            # We need `i` moves on the left branch and `n - i` on the right branch
            # to make sure the corresponding subtrees don't overlap and everything
            # is covered exactly once
            thiscost = sum([
                cost(n - i, leftmost),    # We still need to finish off our walk to the left
                i * left_cost[leftmost],  # The cost for the moves on the left branch
                cost(i, False),           # The tree on the right side, now definitely not leftmost
                right_cost[leftmost] + (n - i - 1) * right_cost[False],  # The cost of moving right, maybe one at the first right cost
            ])
            # If a new lower cost has been found, update values
            if thiscost < c:
                c = thiscost
                checkpoints[leftmost][n] = i
        return c

    def convert(n, checkpoints):
        """
        Given a list of checkpoints, convert this to a list of
        the number of doublings to compute and keep before 
        pushing everything through an isogeny. This forces the
        output to match the more usual implementation, e.g.
        https://crypto.stackexchange.com/a/58377

        Warning! Everything about this function is very hacky, but does the job!
        """
        kernels = [n]
        doubles = []
        leftmost = 1

        # We always select the last point in our kernel
        while kernels != []:
            point = kernels[-1]
            if point == 1:
                # Remove this point and push everything through the isogeny
                kernels.pop()
                kernels = [k - 1 for k in kernels]
                leftmost = 0
            else:
                # checkpoints tells us to double this d times
                d = checkpoints[leftmost][point]
                # Remember that we did this
                doubles.append(d)
                kernels.append(point - d)
        return doubles

    # Compute the cost and populate the checkpoints
    c = cost(n, True)

    # Use the checkpoints to compute the list
    l = convert(n, checkpoints)

    return l


def optimised_strategy(n):
    """"""
    M, S, I = COST[n]['M'], COST[n]['S'], COST[n]['I']
    # Define the costs and initalise the nodes which we store during doubling
    pre_cost = (4*S + 21*M + 1*I, 4*S + 12*M)                                     # (precious_pre_cost, new_pre_cost)
    dbl_cost = ((8*S + 6*M, 12*S + 12*M), (8*S + 8*M, (8*S + 14*M, 6*S + 22*M)))  # (previous_dbl_cost, new_dbl_cost)
    img_cost = ((4*S + 3*M, 18*S + 82*M + 1*I), (4*S + 4*M, 18*S + 81*M))         # (previous_img_cost, new_img_cost)
    cod_cost = ((8*S + 23*M + 1*I, 8*S + 13*M + 1*I), (8*S + 9*M, 8*S + 4*M))     # (previous_cod_cost, new_cod_cost)
    checkpoints = ({}, {})  # (inner, left edge)

    def PREcost(flag):
        return pre_cost[flag]
    
    def DBLcost(n, flag, leftmost):
        cost_per_dbl = dbl_cost[flag][leftmost]
        if flag and leftmost:
            cost = n * cost_per_dbl[0] + cost_per_dbl[1]
        else:
            cost = n * cost_per_dbl
        return cost
    
    def EVALcost(n, Lflag, leftmost):
        cost = 0
        for i in range(n):
            cost = cost + img_cost[Lflag[i]][leftmost]
        return cost
    
    def CODOMAINcost(flag, leftmost, precomp):
        if flag:
            cost = cod_cost[flag][leftmost]
        else:
            cost = cod_cost[flag][precomp]
        return cost

    @functools.cache
    def cost(n, flag, leftmost, precomp):
        """
        The minimal cost to get to all children of a height `n` tree.
        If `flag` is true, we use the inverse-eliminated method.
        If `leftmost` is true, we're still on the leftmost edge of the "outermost" tree.
        If `precomp` is false, we need to consider the cost of precomputing.
        Specially, if `leftmost` is true, we enforce `precomp` to be true.

        To get the mincost, we need compare the result of `cost(n, False, True, True)` and `cost(n, True, True, True)`

        Updates a global "Check points" which are the points along a branch which we 
        keep for later
        """
        if n <= 1:
            return CODOMAINcost(flag, leftmost, precomp), (flag, )  # cost of codomain computing and flag

        mincost = float("inf")
        for i in range(1, n):  # where to branch off
            thiscost = 0
            if not precomp:
                thiscost = PREcost(flag)  # cost of precomputing
            # We need `i` moves on the left branch and `n - i` on the right branch
            # to make sure the corresponding subtrees don't overlap and everything
            # is covered exactly once
            thiscost = thiscost + 2 * DBLcost(i, flag, leftmost)              # cost of doubling, need to double 2 points
            Lcost, Lflag = cost(n-i, flag, leftmost, True)                    # cost of the left branch
            thiscost = thiscost + Lcost + 2 * EVALcost(n-i, Lflag, leftmost)  # cost of evaluating images
            RcostOLD, RflagOLD = cost(i, False, False, False)                 # cost of the right branch using previous method
            RcostNEW, RflagNEW = cost(i, True, False, False)                  # cost of the right branch using new method
            if RcostOLD < RcostNEW:
                Rcost, Rflag = RcostOLD, RflagOLD
            else:
                Rcost, Rflag = RcostNEW, RflagNEW
            thiscost = thiscost + Rcost
            thisflag = Lflag + Rflag
            
            if thiscost < mincost:
                mincost = thiscost
                minflag = thisflag
                checkpoints[leftmost][n] = i
            
        return mincost, minflag

    def convert(n, checkpoints):
        """
        Given a list of checkpoints, convert this to a list of
        the number of doublings to compute and keep before 
        pushing everything through an isogeny. This forces the
        output to match the more usual implementation, e.g.
        https://crypto.stackexchange.com/a/58377

        Warning! Everything about this function is very hacky, but does the job!
        """
        kernels = [n]
        doubles = []
        leftmost = 1

        # We always select the last point in our kernel
        while kernels != []:
            point = kernels[-1]
            if point == 1:
                # Remove this point and push everything through the isogeny
                kernels.pop()
                kernels = [k - 1 for k in kernels]
                leftmost = 0
            else:
                # checkpoints tells us to double this d times
                d = checkpoints[leftmost][point]
                # Remember that we did this
                doubles.append(d)
                kernels.append(point - d)
        return doubles

    # Compute the cost and populate the checkpoints
    costOLD, flagOLD = cost(n, False, True, True)
    checkpointsOLD = tuple(d.copy() for d in checkpoints)
    costNEW, flagNEW = cost(n, True, True, True)
    if costOLD < costNEW:
        mincost = costOLD
        flag = flagOLD
        checkpoints = checkpointsOLD
    else:
        mincost = costNEW
        flag = flagNEW

    # Use the checkpoints to compute the list
    doubles = convert(n, checkpoints)

    return {"doubles": doubles, "flag": flag, }
# fmt: on
# n = length of chain
# M = Cost of Multiplication
# S = Cost of Squaring
# I = Cost of Inversion
# data = [n, M, S, I]

# NOTE:
# Costings are computed from the msi.rs benchmark
# where the cost is the time in ns for one operation
data_sml = [126, 74, 52, 3314]
data_med = [208, 188, 153, 5939]
data_big = [632, 2717, 2265, 54_823]

COST = {
    126: {'M': 74, 'S': 52, 'I': 3314},
    124: {'M': 74, 'S': 52, 'I': 3314},
    208: {'M': 188, 'S': 153, 'I': 5939},
    206: {'M': 188, 'S': 153, 'I': 5939},
    632: {'M': 2717, 'S': 2265, 'I': 54823},
    630: {'M': 2717, 'S': 2265, 'I': 54823},
}

strat_sml = optimised_strategy(data_sml[0])
print(strat_sml)
strat_med = optimised_strategy(data_med[0])
print(strat_med)
strat_big = optimised_strategy(data_big[0])
print(strat_big)
