# Sage has a built-in rref() method for matrices, but sometimes when
# teaching it's nice to get a list of the row operations you would do to
# put the matrix into rref -- or, equivalent, a list of the
# corresponding elementary matrices.
#
# This script contains two functions along those lines:
#
#    rref(): returns a list of elementary matrices that put the matrix into rref
#    rref_steps(): returns a list of strings that describe what row ops to do
#
# The others are just helper functions.
#
# Enjoy.
#
# -Dan Drake, '@'.join(['ddrake', 'member.ams.org'])

def swap_rows(i, j, n):
    ret = identity_matrix(n)
    ret[i, i] = 0
    ret[j, j] = 0
    ret[i, j] = 1
    ret[j, i] = 1
    return ret

def add_c_times_row_i_to_row_j(c, i, j, n):
    ret = identity_matrix(parent(c), n)
    ret[j, i] = c
    return ret

def scale_row_i_by_c(i, c, n):
    ret = identity_matrix(parent(c), n)
    ret[i, i] = c
    return ret

def identify_elem(E):
    """
    Returns a string describing what the elementary matrix E does.

    String is TeX-friendly and indices are 1-based.
    """
    for i in range(E.nrows()):
        for j in range(i + 1, E.nrows()):
            if E[i, j] != 0:
                if E[i, i] == 0:
                    return 'swap rows {0} and {1}'.format(i+1, j+1)
                else:
                    return 'add ${0}$ times row {1} to row {2}'.format(E[i, j], j+1, i+1)

    # no entries above diagonal...
    for i in range(E.nrows()):
        for j in range(i):
            if E[i, j] != 0:
                return 'add ${0}$ times row {1} to row {2}'.format(E[i, j], j+1, i+1)

    # E is diagonal; what row does it scale?
    for i in range(E.nrows()):
        if E[i, i] != 1:
            return 'scale row {0} by ${1}$'.format(i+1, E[i, i])

    # if we get here, we got the identity
    raise ValueError, "the identity matrix is probably a valid elementary matrix, but it's certainly not very useful."

def rref_steps(A):
    """
    Return a list describing the steps one would take to reduce `A` to its rref.
    """
    return [identify_elem(e) for e in reversed(rref(A))]

def rref(A, start_row=0, start_col=0):
    """
    Return `ems`, a list of elementary matrices that turn A into its
    rref. Assumes an exact ring (actually, it assumes QQ), and is not
    optimized; this is intended to represent the basic algorithm
    described in textbooks.

    The returned list satisfies

    A.rref() = prod(reversed(ems)) * A
    """
    ems = []
    n = A.nrows()
    pivot_row = -1
    pivot_col = -1

    # find the pivot: find first nonzero column, and nonzero element in
    # that row.
    for j, col in list(enumerate(A.columns()))[start_col:]:
        try:
            pivot_row = [i for i, e in enumerate(col) if e != 0 and i >= start_row][0]
            pivot_col = j
            break
        except IndexError:
            pass
    
    if pivot_col == -1:
        # zero matrix, it's in rref already.
        return []

    # swap rows so the pivot is at the top.
    if pivot_row > start_row:
        E = swap_rows(start_row, pivot_row, n)
        ems.insert(0, E)
        A = E*A

    # scale pivot to 1
    if A[start_row, pivot_col] != 1:
        E = scale_row_i_by_c(start_row, 1/A[start_row, pivot_col], n)
        ems.insert(0, E)
        A = E*A

    # use pivot to clear its column.
    for row, entry in [(i, e) for i, e in enumerate(A.column(pivot_col))
                       if e != 0 and i != start_row]:
        E = add_c_times_row_i_to_row_j(-entry, start_row, row, n)
        ems.insert(0, E)
        A = E*A

    # recurse
    return rref(A, start_row + 1, pivot_col + 1) + ems

def test(n=10):
    for _ in range(n):
        nrows = randrange(1,15)
        ncols = randrange(1,15)
        A = random_matrix(QQ, nrows, ncols)
        ems = rref(A)
        if prod(ems) * A != A.rref():
            print 'UH OH:'
            print A
            return False
    return True

def row_steps(A):
    ems   = rref(A)
    steps = rref_steps(A)
    
    B = A

    ems.reverse()

    print "Input matrix:"
    print A
    print '\n'

    for i in range(len(ems)):
        B = ems[i]*B
        print steps[i]
        print B
        print '\n'
    return None
