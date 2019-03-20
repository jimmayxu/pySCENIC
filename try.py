from multiprocessing import Pool, cpu_count





import os
import glob
import pickle
import pandas as pd
import numpy as np

from dask.diagnostics import ProgressBar

from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2

from pyscenic.rnkdb import FeatherRankingDatabase as RankingDatabase
from pyscenic.utils import modules_from_adjacencies, load_motifs
from pyscenic.prune import prune2df, df2regulons
from pyscenic.aucell import aucell

import seaborn as sns

DATA_FOLDER="/home/ubuntu/PyCharm/data"
RESOURCES_FOLDER="/home/ubuntu/PyCharm/resources"
DATABASE_FOLDER = "/home/ubuntu/PyCharm/databases/"
SCHEDULER="123.122.8.24:8786"

DATABASES_GLOB = os.path.join(DATABASE_FOLDER, "mm9-*.mc9nr.feather")
MOTIF_ANNOTATIONS_FNAME = os.path.join(RESOURCES_FOLDER, "motifs-v9-nr.mgi-m0.001-o0.0.tbl")

MM_TFS_FNAME = os.path.join(RESOURCES_FOLDER, 'mm_tfs.txt')
SC_EXP_FNAME = os.path.join(RESOURCES_FOLDER, "GSE60361_C1-3005-Expression.txt")

ADJACENCIES_FNAME = os.path.join(DATA_FOLDER, "adjacencies.tsv")
MODULES_FNAME = os.path.join(DATA_FOLDER, "modules.p")
MOTIFS_FNAME = os.path.join(DATA_FOLDER, "motifs.csv")
REGULONS_FNAME = os.path.join(DATA_FOLDER, "regulons.p")

N_SAMPLES = 500


ex_matrix = pd.read_csv(SC_EXP_FNAME, sep='\t', header=0, index_col=0).T
ex_matrix.shape

tf_names = load_tf_names(MM_TFS_FNAME)

db_fnames = glob.glob(DATABASES_GLOB)
def name(fname):
    return os.path.basename(fname).split(".")[0]
dbs = [RankingDatabase(fname=fname, name=name(fname)) for fname in db_fnames]
dbs

# Phase I: Inference of co-expression modules
if __name__ == '__main__':
    adjacencies = grnboost2(ex_matrix, tf_names=tf_names, verbose=True)

adjacencies.head()
adjacencies.to_csv(ADJACENCIES_FNAME, index=False, sep='\t')


modules = list(modules_from_adjacencies(adjacencies, ex_matrix))
with open(MODULES_FNAME, 'wb') as f:
    pickle.dump(modules, f)

#with open(MODULES_FNAME, 'rb') as f:
#    modules = pickle.load(f)



# Phase II: Prune modules for targets with cis regulatory footprints (aka RcisTarget)

df = prune2df(dbs, modules, MOTIF_ANNOTATIONS_FNAME)



with open(REGULONS_FNAME, 'wb') as f:
    pickle.dump(regulons, f)

#with open(REGULONS_FNAME, 'rb') as f:
#    regulons = pickle.load(f)




# Phase III: Cellular regulon enrichment matrix (aka AUCell)

auc_mtx = aucell(ex_matrix, regulons, num_workers=1)





