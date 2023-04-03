import datetime
import hashlib
import json
import os
import os.path
import pandas as pd
from shutil import copytree
from shutil import rmtree
import pathlib
import json
import yaml
import hubmapbags
from pathlib import Path
import pandas as pd
from tqdm import tqdm
from warnings import warn as warning
from datetime import datetime
import shutil
import warnings
from pathlib import Path
import hubmapinventory
import magic  # pyton-magic
import gzip
import numpy as np
import pandas as pd
import tabulate
from tqdm import tqdm

###############################################################################################################
#DISCLAIMER: @icaoberg this code is super alpha. Please be kind.
try:
    from pandas.core.common import SettingWithCopyWarning

    warnings.simplefilter(action="ignore", category=SettingWithCopyWarning)
except:
    warnings.filterwarnings("ignore")


###############################################################################################################
def __pprint(msg: str):
    row = len(msg)
    h = "".join(["+"] + ["-" * row] + ["+"])
    result = "\n" + h + "\n" "|" + msg + "|" "\n" + h
    print(result)


def __update_dataframe(
    dataset: pd.DataFrame, temp: pd.DataFrame, key: str
) -> pd.DataFrame:
    for index, datum in temp.iterrows():
        dataset.loc[index, key] = temp.loc[index, key]
    return dataset


###############################################################################################################
def submission(
    hubmap_ids: list[str],
    dbgap_study_id: str | None,
    token: str | None,
    ncores: int = 2,
    compute_uuids: bool = False,
    debug: bool = False,
) -> bool:
    """
    Main function that creates a dbGaP submission
    """

    for hubmap_id in hubmap_ids:
        hubmapinventory.inventory.create(hubmap_id = hubmap_id,
            dbgap_study_id = dbgap_study_id,
            token = token,
            ncores = ncores,
            compute_uuids = compute_uuids,
            recompute_file_extension = False,
            debug = debug,
        )

###############################################################################################################
#DISCLAIMER: @icaoberg this code is super alpha. Please be kind.
# remove submission folder if it exists
directory = dbgap_study_id
p = pathlib.Path( directory )
if p.exists() and p.is_dir():
    print( 'Removing existing folder ' + directory )
    rmtree(p)
#result = copytree( 'dbgap-submission-scripts/templates', directory )

###############################################################################################################

