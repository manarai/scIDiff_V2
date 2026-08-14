from .decompose import jacobian_modes, jacobian_modes_svd
from .koopman import koopman_modes


def decompose_archetypes(J_tensor, *, method: str = "snmf", rank: int = 5, **kwargs):
    """
    Dispatch a Jacobian-tensor decomposition to a chosen backend.

    Kept as a light wrapper so a caller who already holds a Jacobian
    tensor can compare backends without going through :func:`fit_drift`.

    Parameters
    ----------
    J_tensor : (T, D, D) or (T, N, D, D) torch.Tensor
        Temporal Jacobian tensor produced by scJDO.
    method   : {"snmf", "koopman", "svd"}
        - ``"snmf"``     — non-negative activations of signed pattern operators.
        - ``"koopman"``  — windowed EDMD in a reduced trajectory basis; returns
          eigenvalue-based growth / oscillation diagnostics as well.
        - ``"svd"``      — plain SVD baseline (legacy).
    rank     : int
        Number of modes / archetypes.
    **kwargs : forwarded to the chosen backend.

    Returns
    -------
    The chosen backend's return value. For ``"koopman"`` with
    ``return_diagnostics=True`` (default here) this is a 4-tuple
    ``(patterns, activations, error, diagnostics)``; the semi-NMF and
    SVD backends return their canonical 3-tuple.
    """
    method = method.lower()
    if method == "snmf":
        return jacobian_modes(J_tensor, rank=rank, **kwargs)
    if method == "koopman":
        kwargs.setdefault("return_diagnostics", True)
        return koopman_modes(J_tensor, rank=rank, **kwargs)
    if method == "svd":
        return jacobian_modes_svd(J_tensor, rank=rank, **kwargs)
    raise ValueError(f"Unknown method {method!r}; choose from 'snmf', 'koopman', 'svd'.")


__all__ = [
    "jacobian_modes",
    "jacobian_modes_svd",
    "koopman_modes",
    "decompose_archetypes",
]
