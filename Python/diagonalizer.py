import numpy as np
from scipy.optimize import linear_sum_assignment
_HAS_HUNGARIAN = True


import numpy as np

try:
    from scipy.optimize import linear_sum_assignment
    _HAS_HUNGARIAN = True
except Exception:
    _HAS_HUNGARIAN = False


def diagonalize_ordered(H, prev_vecs=None):
    """
    Diagonalize a Hermitian matrix and ensure consistent ordering of eigenstates.

    Parameters
    ----------
    H : (n, n) ndarray
        Hermitian matrix (Hamiltonian).
    prev_vecs : (n, n) ndarray, optional
        Reference eigenvectors from the previous step. Columns are eigenvectors.
        If provided, the returned eigenvectors are reordered and rephased
        to maximize overlap with prev_vecs.

    Returns
    -------
    eigvals : (n,) ndarray
        Eigenvalues in the same order as matched eigenvectors.
    eigvecs : (n, n) ndarray
        Columns are eigenvectors ordered and phase-aligned with prev_vecs.

    Notes
    -----
    - Uses numpy.linalg.eigh for Hermitian matrices.
    - If prev_vecs is None, eigenvectors are phase-fixed deterministically.
    - Requires SciPy for optimal matching; falls back to greedy matching otherwise.
    """
    H = np.array(H, dtype=np.complex128)
    n = H.shape[0]
    if H.shape[1] != n:
        raise ValueError("Matrix must be square")

    eigvals, eigvecs = np.linalg.eigh(H)

    def _fix_phases(vecs):
        """Make each column's largest component real and positive."""
        vecs = np.array(vecs, dtype=np.complex128)
        for j in range(vecs.shape[1]):
            col = vecs[:, j]
            k = np.argmax(np.abs(col))
            if col[k] != 0:
                phase = col[k] / abs(col[k])
                vecs[:, j] *= np.conj(phase)
        return vecs

    if prev_vecs is None:
        return eigvals, _fix_phases(eigvecs)

    prev_vecs = np.array(prev_vecs, dtype=np.complex128)
    if prev_vecs.shape != (n, n):
        raise ValueError("prev_vecs must have shape (n, n)")

    # Overlap matrix
    overlap = np.abs(prev_vecs.conj().T @ eigvecs)

    # Determine permutation maximizing overlaps
    if _HAS_HUNGARIAN:
        cost = -overlap
        row_ind, col_ind = linear_sum_assignment(cost)
        perm = np.empty(n, dtype=int)
        perm[row_ind] = col_ind
    else:
        perm = np.full(n, -1, dtype=int)
        used_cols = np.zeros(n, bool)
        for _ in range(n):
            i, j = divmod(np.argmax(overlap), n)
            if overlap[i, j] < 0:
                break
            if not used_cols[j]:
                perm[i] = j
                used_cols[j] = True
            overlap[i, :] = -1
            overlap[:, j] = -1
        if np.any(perm < 0):
            free_cols = np.where(~used_cols)[0]
            perm[perm < 0] = free_cols[:np.sum(perm < 0)]

    eigvals = eigvals[perm]
    eigvecs = eigvecs[:, perm]

    # Fix phase relative to previous vectors
    for i in range(n):
        ov = np.vdot(prev_vecs[:, i], eigvecs[:, i])
        if ov != 0:
            eigvecs[:, i] *= np.conj(ov / abs(ov))
        else:
            eigvecs[:, i] = _fix_phases(eigvecs[:, i:i+1])[:, 0]

    return eigvals, eigvecs


if __name__ == "__main__":

    matrix = np.array([[4.0, -30.0, 60.0, -35.0],
                    [-30.0, 300.0, -675.0, 420.0],
                    [60.0, -675.0, 1620.0, -1050.0],
                    [-35.0, 420.0, -1050.0, 700.0]])

    energies, psi = diagonalize_ordered(matrix)

    print(energies, "\n", psi)
