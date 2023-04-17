import pathlib
import warnings
from pathlib import Path
from shutil import rmtree

import hubmapbags
import hubmapinventory
import magic  # pyton-magic
import numpy as np
import pandas as pd
import requests
import tabulate
import yaml
from tqdm import tqdm

###############################################################################################################
# DISCLAIMER: @icaoberg this code is super alpha. Please be kind.
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
    dbgap_study_id: str,
    token: str,
) -> bool:
    """
    Main function that creates a dbGaP submission
    """

    directory = dbgap_study_id
    p = pathlib.Path(directory)
    if p.exists() and p.is_dir():
        print("Removing existing folder " + directory)
        rmtree(p)

    Path(directory).mkdir()

    columns = [
        "sample_id",
        "donor_uuid",
        "donor_hubmap_id",
        "direct_sample_uuid",
        "direct_sample_hubmap_id",
        "organ_uuid",
        "organ_hubmap_id",
        "organ_type",
        "direct_sample_type",
        "dataset_metadata",
        "donor_metadata",
    ]

    df = pd.DataFrame(columns=columns)
    df["sample_id"] = hubmap_ids

    for index, dataset in tqdm(df.iterrows()):
        pmetadata = hubmapbags.apis.get_provenance_info(
            dataset["sample_id"], instance='prod', token=token
        )

        try:
            df.loc[index, "donor_uuid"] = pmetadata["donor_uuid"][0]
        except Exception as e:
            print(e)
            print(pmetadata["donor_uuid"])

        try:
            df.loc[index, "donor_hubmap_id"] = pmetadata["donor_hubmap_id"][0]
        except Exception as e:
            print(e)
            print(pmetadata["donor_hubmap_id"])

        df.loc[index, "direct_sample_uuid"] = pmetadata["first_sample_uuid"][0]
        df.loc[index, "direct_sample_type"] = pmetadata["first_sample_type"][0]
        df.loc[index, "direct_sample_hubmap_id"] = pmetadata["first_sample_hubmap_id"][
            0
        ]

        try:
            df.loc[index, "organ_uuid"] = pmetadata["organ_uuid"][0]
        except Exception as e:
            print(e)
            print(pmetadata["organ_uuid"])

        try:
            df.loc[index, "organ_hubmap_id"] = pmetadata["organ_hubmap_id"][0]
        except Exception as e:
            print(e)
            print(pmetadata["organ_hubmap_id"])

        try:
            df.loc[index, "organ_type"] = pmetadata["organ_type"][0]
        except Exception as e:
            print(e)
            print(pmetadata["organ_type"])

        metadata = hubmapbags.apis.get_dataset_info(
            dataset["sample_id"], instance='prod', token=token
        )

        try:
            df.loc[index, "donor_uuid"] = pmetadata.get("donor_uuid")[0]
        except Exception as e:
            print(e)
            print(pmetadata.get("donor_uuid"))

        try:
            df.loc[index, "donor_hubmap_id"] = pmetadata.get("donor_hubmap_id")[0]
        except Exception as e:
            print(e)
            print(pmetadata.get("donor_hubmap_id"))

    __create_donor_metadata(df, token, directory)
    __create_sample_attributes(df, token, directory)
    __create_sample_mapping(df, token, directory)
    __get_spreadhsheets(directory)

    return True


###############################################################################################################
def __create_donor_metadata(df: pd.DataFrame, token: str, directory: str) -> None:
    donor = df[["donor_hubmap_id", "donor_uuid"]]
    donor = donor.drop_duplicates(subset=["donor_hubmap_id"])

    donor["sex"] = None
    for index, datum in tqdm(donor.iterrows()):
        metadata = hubmapbags.apis.get_entity_info(
            datum["donor_hubmap_id"], token=token, instance="prod"
        )
        if "living_donor_data" in metadata["metadata"].keys():
            for info in metadata["metadata"]["living_donor_data"]:
                if info["grouping_concept_preferred_term"] == "Sex":
                    donor.loc[index, "sex"] = info["preferred_term"]
        else:
            for info in metadata["metadata"]["organ_donor_data"]:
                if info["grouping_concept_preferred_term"] == "Sex":
                    donor.loc[index, "sex"] = info["preferred_term"]

        if donor.loc[index, "sex"] == "Male":
            donor.loc[index, "sex"] = 1
        else:
            donor.loc[index, "sex"] = 2

        donor.loc[index, "subject_source"] = "HuBMAP"

    donor = donor.drop("donor_uuid", axis=1)
    donor["SOURCE_SUBJECT_ID"] = donor["donor_hubmap_id"]
    donor["consent"] = 1
    donor = donor.rename(
        columns={
            "donor_hubmap_id": "SUBJECT_ID",
            "consent": "CONSENT",
            "sex": "SEX",
            "subject_source": "SUBJECT_SOURCE",
        }
    )
    donor = donor.reindex(
        columns=["SUBJECT_ID", "CONSENT", "SEX", "SUBJECT_SOURCE", "SOURCE_SUBJECT_ID"]
    )
    donor.to_csv(f"{directory}/2a_SubjectConsent_DS.txt", index=False, sep="\t")

###############################################################################################################
def __create_sample_attributes(df: pd.DataFrame, token: str, directory: str):
    URL = "https://raw.githubusercontent.com/hubmapconsortium/search-api/main/src/search-schema/data/definitions/enums/organ_types.yaml"

    temp_file = Path("/tmp/organ_types.yaml")
    if temp_file.exists():
        temp_file.unlink()

    response = requests.get(URL)
    temp_file.write_bytes(response.content)

    with open("/tmp/organ_types.yaml") as file:
        organ_types = yaml.load(file, Loader=yaml.FullLoader)

    sample_attributes = df[["sample_id"]]
    analyte_class = []

    sample_attributes["BODY_SITE"] = None
    for index, datum in tqdm(sample_attributes.iterrows()):
        metadata = hubmapbags.apis.get_dataset_info(
            datum["sample_id"], token=token, instance='prod'
        )

        if datum["sample_id"] == "HBM347.RFGL.437":
            analyte_class.append("DNA")
        elif datum["sample_id"] == "HBM773.WCXC.264":
            analyte_class.append("RNA")
        elif "ingest_metadata" in metadata.keys():
            analyte_class.append(
                metadata["ingest_metadata"]["metadata"]["analyte_class"]
            )
        else:
            print(datum["sample_id"])

        sample_attributes.loc[index, "BODY_SITE"] = df.loc[index, "organ_type"]

    sample_attributes["ANALYTE_TYPE"] = analyte_class
    sample_attributes["IS_TUMOR"] = "N"
    sample_attributes = sample_attributes.rename(columns={"sample_id": "SAMPLE_ID"})
    sample_attributes = sample_attributes.reindex(
        columns=["SAMPLE_ID", "BODY_SITE", "ANALYTE_TYPE", "IS_TUMOR"]
    )
    sample_attributes.to_csv(
        directory + "/6a_SampleAttributes_DS.txt", index=False, sep="\t"
    )


###############################################################################################################
def __create_sample_mapping(df: pd.DataFrame, token: str, directory: str):
    sample_mapping = df[["donor_hubmap_id", "sample_id"]]
    sample_mapping = sample_mapping.rename(
        columns={"donor_hubmap_id": "SUBJECT_ID", "sample_id": "SAMPLE_ID"}
    )
    sample_mapping.to_csv(f'{directory}/3a_SSM_DS.txt', index=False, sep="\t")

###############################################################################################################
def __get_spreadhsheets(directory: str):
    filename = '2b_SubjectConsent_DD.xlsx'
    URL = "https://github.com/hubmapconsortium/py-hubmap-dbgap/blob/master/files/2b_SubjectConsent_DD.xlsx?raw=true"
    temp_file = Path(f'{directory}/{filename}')
    response = requests.get(URL)
    temp_file.write_bytes(response.content)

    filename = '3b_SSM_DD.xlsx'
    URL = "https://github.com/hubmapconsortium/py-hubmap-dbgap/blob/master/files/3b_SSM_DD.xlsx?raw=true"
    temp_file = Path(f'{directory}/{filename}')
    response = requests.get(URL)
    temp_file.write_bytes(response.content)

    filename = '6b_SampleAttributes_DD.xlsx'
    URL = "https://github.com/hubmapconsortium/py-hubmap-dbgap/blob/master/files/6b_SampleAttributes_DD.xlsx?raw=true"
    temp_file = Path(f'{directory}/{filename}')
    response = requests.get(URL)
    temp_file.write_bytes(response.content)



