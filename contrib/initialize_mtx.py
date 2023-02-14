from pyliger import create_liger, read_10X_h5, normalize, scale_not_center, select_genes
from pyliger.factorization._online_iNMF import _get_scale_data
from scipy import io, sparse
import numpy as np
from os import listdir

titlelist = []
h5list = []
for file in listdir('../../repos/theAeon/pyliger/h5'):
    titlelist.append(file)
    h5list.append(read_10X_h5(sample_dir='../../repos/theAeon/pyliger/h5',
                  sample_name=file.split(".")[0], file_name=file, backed=True))
for i in range(0, len(h5list)//2):
    j = i + len(h5list)//2
    titlepair = (titlelist[i], titlelist[j])
    mouse_cortex = create_liger([h5list[i], h5list[j]])
    normalize(mouse_cortex)
    select_genes(mouse_cortex, var_thresh=0.2, do_plot=False)
    scale_not_center(mouse_cortex)
    outdata = [_get_scale_data(adata)['scale_data'] for adata in mouse_cortex.adata_list for i in range(2)]
    outdata = [mat.value.toarray() for mat in outdata]
    outdata = [sparse.coo_array(np.hstack((mat, np.zeros(mat.shape))).T) for mat in outdata]
    for i in range(len(titlepair)):
        io.mmwrite( "%s.mtx" % titlepair[i], outdata[i])
    densedata = np.random.uniform(0, 2, (10, 2 * len(mouse_cortex.var_genes)))
    io.mmwrite("%s.dense.mtx" % titlepair[0], densedata)
