from urllib import request
from tqdm.auto import tqdm

class DownloadProgressBar(tqdm):
    def update_to(self, b=1, bsize=1, tsize=None):
        if tsize is not None:
            self.total = tsize
        self.update(b * bsize - self.n)

def install_databases(install=['human', 'mouse']):   
    database_human = "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-tss-centered-10kb-7species.mc9nr.feather"
    database_mouse = "https://resources.aertslab.org/cistarget/databases/mus_musculus/mm9/refseq_r45/mc9nr/gene_based/mm9-500bp-upstream-7species.mc9nr.genes_vs_motifs.rankings.feather"
    if "human" in install:
        with DownloadProgressBar(unit='B', unit_scale=True, miniters=1, desc=database_human.split('/')[-1]) as t:
            request.urlretrieve(database_human, filename="databases/hg19-tss-centered-10kb-7species.mc9nr.feather", reporthook=t.update_to)
    if "mouse" in install:
        with DownloadProgressBar(unit='B', unit_scale=True, miniters=1, desc=database_mouse.split('/')[-1]) as t:
            request.urlretrieve(database_mouse, filename="databases/mm9-500bp-upstream-7species.mc9nr.genes_vs_motifs.rankings.feather", reporthook=t.update_to)