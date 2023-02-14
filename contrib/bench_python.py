import itertools
import pathlib
from os import PathLike
from typing import TypeGuard, Union
import numpy as np
from scipy import io, sparse, optimize as opt
from pyliger.factorization._utilities import nla, nnlsm_blockpivot
from nnlsm_activeset import nnlsm_activeset
import sys

def parse_glob(p: pathlib.Path) -> list[pathlib.Path]:
    if p.is_dir():
        return list(p.glob('*.mtx'))
    else:
        return [p]


def safe_mmread(mt: PathLike) -> Union[sparse.csr_array, None]:
    try:
        return sparse.csr_array(io.mmread(mt).T)
    except ValueError as e:
        raise e


def is_csr_array(k: sparse.csr_array | None) -> TypeGuard[sparse.csr_array]:
    return isinstance(k, sparse.csr_array)


def _test_nnlsm(A_in, B_in):
    ### 1. Initialization (W, V_i, H_i)
    A = A_in.toarray()
    B = B_in.toarray()
    X_org = A.T.dot(B)
    # B = np.random.rand(m,k)
    # A = np.random.rand(m,n/2)
    # A = np.concatenate((A,A),axis=1)
    # A = A + np.random.rand(m,n)*0.01
    # B = np.random.rand(m,k)
    import time
    start = time.time()
    C1, info = nnlsm_blockpivot(A, B)
    elapsed2 = time.time() - start
    rel_norm2 = nla.norm(C1 - X_org) / nla.norm(X_org)
    print('nnlsm_blockpivot:    ', 'OK  ' if info[0] else 'Fail',
          'elapsed:{0:.4f} error:{1:.4e}'.format(elapsed2, rel_norm2))
    start = time.time()
    C2, info = nnlsm_activeset(A, B)
    elapsed1 = time.time() - start
    rel_norm1 = nla.norm(C2 - X_org) / nla.norm(X_org)
    print('nnlsm_activeset:     ', 'OK  ' if info[0] else 'Fail',
          'elapsed:{0:.4f} error:{1:.4e}'.format(elapsed1, rel_norm1))
    start = time.time()
    C3 = np.zeros(X_org.shape)
    for i in iter(range(0, X_org.shape[1])):
        res = opt.nnls(A, B[:, i])
        C3[:, i] = res[0]
    elapsed3 = time.time() - start
    rel_norm3 = nla.norm(C3 - X_org) / nla.norm(X_org)
    print ('scipy.optimize.nnls: ', 'OK  ',\
        'elapsed:{0:.4f} error:{1:.4e}'.format(elapsed3, rel_norm3))


benchlist = []
assert len(sys.argv) == 3
for i in sys.argv[1:]:
    benchlist.append([path.resolve() for path in parse_glob(pathlib.Path(i))])
benchlist = [safe_mmread(m) for m in itertools.chain.from_iterable(benchlist)]
benchlist = list(filter(is_csr_array, benchlist))
_test_nnlsm(benchlist[0], benchlist[1].T)
