from urllib import request
from pathlib import Path
import shutil
from xmlrpc.client import Boolean
from tqdm.auto import tqdm
import os
import glob

class DownloadProgressBar(tqdm):
    global MAIN_PATH
    MAIN_PATH = os.path.dirname(__file__)
    
    def update_to(self, b=1, bsize=1, tsize=None):
        if tsize is not None:
            self.total = tsize
        self.update(b * bsize - self.n)
        

def db_exists(dataset)-> Boolean:
    dataset_human = os.path.join(MAIN_PATH, "databases/human")
    dataset_mouse = os.path.join(MAIN_PATH, "databases/mouse")
    
    if dataset == "human":
        return os.path.exists(dataset_human)
    if dataset == "mouse":
        return os.path.exists(dataset_mouse)


def download_db():
    print("Do you want to install datasets for (human/mouse/both)?")
    answer=input()
    dataset_human = "https://wolken.zbh.uni-hamburg.de/index.php/s/DZ3XT76Z7xYAT64/download/human.zip"
    dataset_mouse = "https://wolken.zbh.uni-hamburg.de/index.php/s/6syqpRTgDLaQs4t/download/mouse.zip"

    if answer == "human" or answer == "both":
        if db_exists("human"):
            print('Human dataset is already downloaded!')
        else:
            filename_ = os.path.join(MAIN_PATH, "databases",os.path.basename(dataset_human))
            with DownloadProgressBar(unit='B', unit_scale=True, miniters=1, desc=dataset_human.split('/')[-1]) as t:
                request.urlretrieve(dataset_human, filename=filename_,
                                    reporthook=t.update_to)

        
    if answer == "mouse" or answer == "both":
        if db_exists("mouse"):
            print('Mouse dataset is already downloaded!')
        else:
            filename_ = os.path.join(MAIN_PATH, "databases",os.path.basename(dataset_mouse))
            with DownloadProgressBar(unit='B', unit_scale=True, miniters=1, desc=dataset_mouse.split('/')[-1]) as t:
                request.urlretrieve(dataset_mouse, 
                                    filename=filename_, 
                                    reporthook=t.update_to)
        
    if answer not in ["human","mouse","both"]:
        print("Invalid answer")
        
#     # extract
    for f in glob.glob(os.path.join(MAIN_PATH,"databases/*.zip")):
        print(f)
        print(os.path.dirname(f))
        out_dir = os.path.dirname(f)
        print(out_dir)
        shutil.unpack_archive(f, out_dir)
        os.remove(f)