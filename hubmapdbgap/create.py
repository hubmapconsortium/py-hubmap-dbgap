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
    instance='prod',
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

    directory = dbgap_study_id
    p = pathlib.Path( directory )
    if p.exists() and p.is_dir():
        print( 'Removing existing folder ' + directory )
        rmtree(p)

    columns = ['sample_id', 'donor_uuid','donor_hubmap_id',\
           'direct_sample_uuid','direct_sample_hubmap_id',\
           'organ_uuid','organ_hubmap_id','organ_type',\
           'direct_sample_type','dataset_metadata','donor_metadata']
    
    df = pd.DataFrame(columns=columns)
    df['sample_id'] = hubmap_ids
        
    for index, dataset in tqdm(df):
        pmetadata = hubmapbags.apis.get_provenance_info( hubmap_id, \
            instance=instance, token=token)
        
        try:
            df.loc[index,'donor_uuid'] = pmetadata['donor_uuid'][0]
        except Exception as e:
            print(e)
            print(pmetadata['donor_uuid'])
            
        try:
            df.loc[index,'donor_hubmap_id'] = pmetadata['donor_hubmap_id'][0]
        except Exception as e:
            print(e)
            print(pmetadata['donor_hubmap_id'])
        
        df.loc[index,'direct_sample_uuid'] = pmetadata['first_sample_uuid'][0]
        df.loc[index,'direct_sample_type'] = pmetadata['first_sample_type'][0]
        df.loc[index,'direct_sample_hubmap_id'] = pmetadata['first_sample_hubmap_id'][0]
        
        try:
            df.loc[index,'organ_uuid'] = pmetadata['organ_uuid'][0]
        except Exception as e:
            print(e)
            print(pmetadata['organ_uuid'])
        
        try:
            df.loc[index,'organ_hubmap_id'] = pmetadata['organ_hubmap_id'][0]
        except Exception as e:
            print(e)
            print(pmetadata['organ_hubmap_id'])
            
        try:
            df.loc[index,'organ_type'] = pmetadata['organ_type'][0]
        except Exception as e:
            print(e)
            print(pmetadata['organ_type'])
        
        metadata = hubmapbags.apis.get_dataset_info( datum['hubmap_id'], instance=instance, token=token )
        
        try:
            df.loc[index,'donor_uuid'] = pmetadata.get('donor_uuid')[0]
        except Exception as e:
            print(e)
            print(pmetadata.get('donor_uuid'))
            
        try:
            df.loc[index,'donor_hubmap_id'] = pmetadata.get('donor_hubmap_id')[0]
        except Exception as e:
            print(e)
            print(pmetadata.get('donor_hubmap_id'))
    
    return df

###############################################################################################################

donor = report[['donor_hubmap_id', 'donor_uuid']]
donor = donor.drop_duplicates(subset=['donor_hubmap_id'])

donor['sex'] = None
for index, datum in tqdm(donor.iterrows()):
    metadata = hubmapbags.apis.get_entity_info( datum['donor_hubmap_id'], token=token, instance='prod' )
    if 'living_donor_data' in metadata['metadata'].keys():
        for info in metadata['metadata']['living_donor_data']:
            if info['grouping_concept_preferred_term'] == 'Sex':
                donor.loc[index,'sex'] = info['preferred_term']
    else:
        for info in metadata['metadata']['organ_donor_data']:
            if info['grouping_concept_preferred_term'] == 'Sex':
                donor.loc[index,'sex'] = info['preferred_term']
                
    if donor.loc[index,'sex'] == 'Male':
        donor.loc[index,'sex'] = 1;
    else:
        donor.loc[index,'sex'] = 2;
        
    donor.loc[index,'subject_source']='HuBMAP'
    
donor = donor.drop('donor_uuid',axis=1)
donor['SOURCE_SUBJECT_ID']=donor['donor_hubmap_id']
donor['consent']=1
donor = donor.rename(columns={'donor_hubmap_id':'SUBJECT_ID','consent':'CONSENT','sex':'SEX', 'subject_source':'SUBJECT_SOURCE'})
donor=donor.reindex(columns=['SUBJECT_ID', 'CONSENT', 'SEX', 'SUBJECT_SOURCE', 'SOURCE_SUBJECT_ID'])
donor.to_csv(directory + '/2a_SubjectConsent_DS.txt', index=False, sep='\t')

donor

###############################################################################################################
with open('search-api/src/search-schema/data/definitions/enums/organ_types.yaml') as file:
    organ_types = yaml.load(file, Loader=yaml.FullLoader)

sample_attributes = report[['hubmap_id']]
analyte_class = []

sample_attributes['BODY_SITE']=None
for index, datum in tqdm(sample_attributes.iterrows()):
    metadata = hubmapbags.apis.get_dataset_info(datum['hubmap_id'], token=token, instance=instance)
    
    if datum['hubmap_id'] == 'HBM347.RFGL.437':
        analyte_class.append('DNA')
    elif datum['hubmap_id'] == 'HBM773.WCXC.264':
        analyte_class.append('RNA')
    elif 'ingest_metadata' in metadata.keys():
        analyte_class.append(metadata['ingest_metadata']['metadata']['analyte_class'])
    else:
        print(datum['hubmap_id'])
    
    sample_attributes.loc[index,'BODY_SITE'] = report.loc[index, 'organ_type']

sample_attributes['ANALYTE_TYPE'] = analyte_class
sample_attributes['IS_TUMOR'] = 'N'
sample_attributes = sample_attributes.rename(columns={'hubmap_id':'SAMPLE_ID'})
sample_attributes=sample_attributes.reindex(columns=['SAMPLE_ID', 'BODY_SITE', 'ANALYTE_TYPE', 'IS_TUMOR'])
sample_attributes.to_csv(directory + '/6a_SampleAttributes_DS.txt', index=False, sep='\t')

###############################################################################################################
def __create_samlpe_mapping(df: pd.DataFrame, directory: str):
    try:
        sample_mapping = report[['donor_hubmap_id','hubmap_id']]
        sample_mapping = sample_mapping.rename(columns={'donor_hubmap_id':'SUBJECT_ID','hubmap_id':'SAMPLE_ID'})
        sample_mapping.to_csv(directory  + '/3a_SSM_DS.txt', index=False, sep='\t')

        return True
    except:
        return False