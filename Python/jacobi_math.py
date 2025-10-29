## module jacobi
"""
lam,x = jacobi(a,tol = 1.0e-9).
    Solution of std. eigenvalue problem [a]{x} = lam{x}
    by Jacobi's method. Returns eigenvalues in vector {lam}
    and the eigenvectors as columns of matrix [x].
"""
from numpy import array, identity, diagonal
from math import sqrt


def jacobi(a, tol=1.0e-9):
    def maxElem(a):  # Find largest off-diag. element a[k,l]
        n = len(a)
        aMax = 0.0
        for i in range(n - 1):
            for j in range(i + 1, n):
                if abs(a[i, j]) >= aMax:
                    aMax = abs(a[i, j])
                    k = i;
                    l = j
        return aMax, k, l

    def rotate(a, p, k, l):  # Rotate to make a[k,l] = 0
        n = len(a)
        aDiff = a[l, l] - a[k, k]
        if abs(a[k, l]) < abs(aDiff) * 1.0e-36:
            t = a[k, l] / aDiff
        else:
            phi = aDiff / (2.0 * a[k, l])
            t = 1.0 / (abs(phi) + sqrt(phi ** 2 + 1.0))
            if phi < 0.0: t = -t
        c = 1.0 / sqrt(t ** 2 + 1.0);
        s = t * c
        tau = s / (1.0 + c)
        temp = a[k, l]
        a[k, l] = 0.0
        a[k, k] = a[k, k] - t * temp
        a[l, l] = a[l, l] + t * temp
        for i in range(k):  # Case of i < k
            temp = a[i, k]
            a[i, k] = temp - s * (a[i, l] + tau * temp)
            a[i, l] = a[i, l] + s * (temp - tau * a[i, l])
        for i in range(k + 1, l):  # Case of k < i < l
            temp = a[k, i]
            a[k, i] = temp - s * (a[i, l] + tau * a[k, i])
            a[i, l] = a[i, l] + s * (temp - tau * a[i, l])
        for i in range(l + 1, n):  # Case of i > l
            temp = a[k, i]
            a[k, i] = temp - s * (a[l, i] + tau * temp)
            a[l, i] = a[l, i] + s * (temp - tau * a[l, i])
        for i in range(n):  # Update transformation matrix
            temp = p[i, k]
            p[i, k] = temp - s * (p[i, l] + tau * p[i, k])
            p[i, l] = p[i, l] + s * (temp - tau * p[i, l])

    n = len(a)
    maxRot = 5 * (n ** 2)  # Set limit on number of rotations
    p = identity(n) * 1.0  # Initialize transformation matrix
    for i in range(maxRot):  # Jacobi rotation loop
        aMax, k, l = maxElem(a)
        if aMax < tol:
            return diagonal(a), p

        rotate(a, p, k, l)

    print('Jacobi method did not converge')
    return None


import numpy as np
import random as rand


def jacobi_diagonalize(sym_arr, offset_threshold=10 ** -8, record=lambda off, cur_diag, cur_steps: None, start=lambda mat: None,
                       stop=lambda res: None
                       ):
    n = len(sym_arr)
    diagonal = sym_arr.astype(float).copy()
    steps = np.eye(n)
    off = offset(diagonal)
    start(diagonal)
    record(off, diagonal, steps)
    while off >= offset_threshold:
        largest_offset_pos = largest_off_index(diagonal)
        diagonalized_small = diagonalize_2by2(small_matr_of(diagonal,
                                                            largest_offset_pos))
        as_n_mat = promote(diagonalized_small[0], largest_offset_pos, n)
        diagonal = np.dot(np.dot(np.transpose(as_n_mat), diagonal), as_n_mat)
        steps = np.dot(steps, as_n_mat)
        off = offset(diagonal)
        record(off, diagonal, steps)
    res = [steps, diagonal, np.transpose(steps)]
    stop(res)
    return res


def generate_matrix(rows, cols, generating_function):
    result = []
    for i in range(rows):
        for j in range(cols):
            result.append(generating_function(i, j))
    return np.array(result, dtype=float).reshape(rows, cols)


def random_diagonal_matrix(size, max_element):
    reflected = {}

    def diag_generator(row, col):
        if row <= col:
            val = rand.random() * max_element * (-1) ** rand.randint(0, 1)
            reflected[(row, col)] = int(val + 1)
            return reflected[(row, col)]
        else:
            return reflected[(col, row)]

    return generate_matrix(size, size, diag_generator)


def largest_off_index(sym_matr):
    max = -1
    pos = [0, 0]
    for i in range(len(sym_matr)):
        for j in range(i, len(sym_matr)):
            off = sym_matr[i][j] ** 2
            if off > max and i != j:
                max = off
                pos = [i, j]
    return pos


def diagonalize_2by2(arr):
    eigs = eig_calc(arr)
    vals = eigs[0]
    vects = eigs[1]
    g = np.array([vects[0], vects[1]])
    d = np.array([[vals[0], 0], [0, vals[1]]])
    return [np.transpose(g), d, g]


def eig_calc(arr):
    a = arr[0][0]
    b = arr[0][1]
    c = arr[1][1]
    eigvals = calc_eigvals(a, b, c)
    eigvects = calc_eigvects(a, b, c, eigvals)
    return [eigvals, eigvects]


def calc_eigvals(a, b, c):
    if b == 0:  ##it is already diagonalized!
        return [a, c]
    eigvals = []
    base = (a + c) / 2
    radius = np.sqrt(a ** 2 - 2 * a * c + 4 * b ** 2 + c ** 2) / 2
    eigvals.append(base + radius)
    eigvals.append(base - radius)
    return eigvals


def calc_eigvects(a, b, c, eigvals):
    if b == 0:  ##already diagonal...return I
        if eigvals[0] == a:
            return [np.array([1.0, 0.0]), np.array([0.0, 1.0])]
        elif eigvals[0] == b:
            return [np.array([0.0, 1.0]), np.array([1.0, 0.0])]
        else:
            raise np.err

    def eigvect_of(l):
        v1 = 1
        v2 = (l - a) / b
        mag = np.sqrt(v1 ** 2 + v2 ** 2)
        return np.array([v1 / mag, v2 / mag])

    return list(map(eigvect_of, eigvals))


def small_matr_of(large, pos):
    a = pos[0]
    b = pos[1]
    return np.array([large[a][a], large[a][b], large[b][a],
                     large[b][b]]).reshape(2, 2)


class OffFinder:
    def __init__(self):
        self.pos = [0, 0]

    def largest_off_index(self, sym_matr):
        pos = self.pos
        if pos[1] >= len(sym_matr) - 1:
            if pos[0] >= len(sym_matr) - 2:
                pos[0] = 0
            else:
                pos[0] = pos[0] + 1
            pos[1] = pos[0] + 1  # Not on diagonal
        else:
            pos[1] = pos[1] + 1
        return pos


def promote(two_by_two, pos, up_dim):
    result = np.eye(up_dim)
    a = pos[0]
    b = pos[1]
    result[a][a] = two_by_two[0][0]
    result[a][b] = two_by_two[0][1]
    result[b][a] = two_by_two[1][0]
    result[b][b] = two_by_two[1][1]
    return result


def offset(matr):
    res = 0
    for i in range(len(matr)):
        for j in range(len(matr)):
            if i != j:
                res += abs(matr[i][j])
    return res


# Testing jabobi module

if __name__ == "__main__":
    ############ Test of jacobi module ##############
    from numpy import array
    from numpy.linalg import eig

    matrix = array([[4.0, -30.0, 60.0, -35.0],
                    [-30.0, 300.0, -675.0, 420.0],
                    [60.0, -675.0, 1620.0, -1050.0],
                    [-35.0, 420.0, -1050.0, 700.0]])

    matrix2 = array([[40.0, 0, 0, 0],
                     [4000.0, 30.0, 0, 0],
                     [210.0, 250.0, 20.0, 0],
                     [30.0, 500.0, 30.0, 10.0]])

    matrix3 = matrix2.transpose()

    eigenvalues, diagonalized_matrix = jacobi(matrix3)
    print("Eigenvalues jacobi:")
    print(eigenvalues)

    eigenvalues_np, eigenvectors_np = eig(matrix3)
    print("Eigenvalues numpy:")
    print(eigenvalues_np)


