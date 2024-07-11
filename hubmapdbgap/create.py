import pathlib
import warnings
from pathlib import Path
from shutil import rmtree
import zipfile
import os
from datetime import datetime
import hubmapbags
import hubmapinventory
import magic  # pyton-magic
import numpy as np
import pandas as pd
import requests
import tabulate
import yaml
from tqdm import tqdm

##########################################################################
# DISCLAIMER: @icaoberg this code is super alpha. Please be kind.
try:
    from pandas.core.common import SettingWithCopyWarning

    warnings.simplefilter(action="ignore", category=SettingWithCopyWarning)
except BaseException:
    warnings.filterwarnings("ignore")


##########################################################################
def __compress_folder(folder_path, output_path):
    with zipfile.ZipFile(output_path, "w", zipfile.ZIP_DEFLATED) as zipf:
        for root, dirs, files in os.walk(folder_path):
            for file in files:
                file_path = os.path.join(root, file)
                arcname = os.path.relpath(file_path, folder_path)
                zipf.write(file_path, arcname)


##########################################################################
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


##########################################################################
def __print_to_file(output_filename, string):
    with open(output_filename, "a") as file:
        file.write(string)
    file.close()


def submission(
    hubmap_ids: list[str],
    dbgap_study_id: str,
    token: str,
    prepend_sample_id: bool,
) -> bool:
    """
    Main function that creates a dbGaP submission
    """

    directory = dbgap_study_id
    p = pathlib.Path(directory)
    if p.exists() and p.is_dir():
        print(f'Removing existing folder "{directory}')
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
    print(f"Adding {str(len(hubmap_ids))} to the main dataframe")
    df["sample_id"] = hubmap_ids

    print("Gathering dataset metadata")
    for index, dataset in df.iterrows():
        print(f'Processing dataset {dataset["sample_id"]}')
        pmetadata = hubmapbags.apis.get_provenance_info(
            dataset["sample_id"], instance="prod", token=token
        )

        try:
            df.loc[index, "donor_uuid"] = pmetadata["donor_uuid"][0]
        except Exception as e:
            print(pmetadata)

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
            dataset["sample_id"], instance="prod", token=token
        )

        dataset_directory = f'/hive/hubmap/data/{metadata["local_directory_rel_path"]}'

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

    print("Gathering donor metadata")
    __create_donor_metadata(df, token, directory)
    print("Gathering sample attributes")
    __create_sample_attributes(df, token, directory)
    print("Creating sample mapping")
    __create_sample_mapping(df, token, directory)
    print("Downloading spreadsheets")
    __get_spreadhsheets(directory)

    print("Processing dataset inventories")
    data = []
    for hubmap_id in tqdm(hubmap_ids):
        metadata = hubmapbags.apis.get_dataset_info(
            hubmap_id, instance="prod", token=token
        )

        ometadata = hubmapbags.apis.get_provenance_info(
            hubmap_id, instance="prod", token=token
        )

        # THE METADATA
        try:
            library_id = (
                f'{metadata["ingest_metadata"]["metadata"]["library_id"]}-{hubmap_id}'
            )
        except BaseException:
            library_id = f"lib-{hubmap_id}"

        title = f'{metadata["dataset_type"][0]} of {ometadata["organ_type"][0]}'

        # library_strategy
        library_strategy = {
            "SNARE-ATACseq2": "ATAC-seq",
            "SNARE-RNAseq2": "RNA-Seq",
            "scRNA-Seq-10x": "RNA-Seq",
            "sciRNAseq": "RNA-Seq",
            "sciATACseq": "ATAC-seq",
            "ATACseq-bulk": "ATAC-seq",
            "scRNA-Seq-10x": "RNA-Seq",
            "WGS": "WGS",
            "bulk-RNA": "RNA-Seq",
            "scRNAseq-10xGenomics-v3": "RNA-Seq",
            "snATACseq": "ATAC-seq",
            "Slide-seq": "RNA-Seq",
            "snRNAseq": "RNA-Seq",
            "snRNAseq-10xGenomics-v3": "RNA-Seq",
            "scRNA-Seq-10x": "RNA-Seq",
            "MUSIC": "OTHER",
        }

        analyte_class = {
            "RNA": "TRANSCRIPTOMIC",
            "DNA": "GENOMIC",
            "DNA + RNA": "OTHER",
        }

        library_source = analyte_class[
            metadata["ingest_metadata"]["metadata"]["analyte_class"]
        ]

        if (
            metadata["dataset_type"][0] == "SNAREseq"
            and metadata["ingest_metadata"]["metadata"]["analyte_class"] == "RNA"
        ):
            library_strategy = "RNA-Seq"
        elif (
            metadata["dataset_type"][0] == "SNAREseq"
            and metadata["ingest_metadata"]["metadata"]["analyte_class"] == "DNA"
        ):
            library_strategy = "ATAC-seq"
        elif (
            metadata["dataset_type"][0] == "sciRNAseq"
            and metadata["ingest_metadata"]["metadata"]["analyte_class"] == "RNA"
        ):
            library_strategy = "RNA-Seq"
        elif (
            metadata["dataset_type"][0] == "sciATACseq"
            and metadata["ingest_metadata"]["metadata"]["analyte_class"] == "RNA"
        ):
            library_strategy = "RNA-Seq"
        elif (
            metadata["dataset_type"][0] == "sciATACseq"
            and metadata["ingest_metadata"]["metadata"]["analyte_class"] == "DNA"
        ):
            library_strategy = "ATAC-seq"
        else:
            # library_strategy = library_strategy[metadata["dataset_type"][0]]
            library_strategy = library_strategy[metadata["dataset_type"]]

        library_layout = {"paired-end": "paired", "paired end": "paired"}
        library_layout = library_layout[
            metadata["ingest_metadata"]["metadata"]["library_layout"]
        ]

        # library_selection
        library_selection = "other"

        # platform
        platform = {"Illumina": "ILLUMINA"}
        platform = platform[
            metadata["ingest_metadata"]["metadata"]["acquisition_instrument_vendor"]
        ]

        # instrument_model
        instrument_model = {
            "NovaSeq": "Illumina NovaSeq 6000",
            "NovaSeq6000": "Illumina NovaSeq 6000",
            "Novaseq6000": "Illumina NovaSeq 6000",
            "NovaSeq 6000": "Illumina NovaSeq 6000",
            "Novaseq 6000": "Illumina NovaSeq 6000",
            "HiSeq": "Illumina HiSeq 4000",
            "HiSeq 4000": "Illumina HiSeq 4000",
            "Nextseq2000": "NextSeq 2000",
            "Novaseq6020": "Illumina NovaSeq 6000",
            "Novaseq6019": "Illumina NovaSeq 6000",
            "Novaseq6018": "Illumina NovaSeq 6000",
            "Novaseq6016": "Illumina NovaSeq 6000",
            "Novaseq6015": "Illumina NovaSeq 6000",
            "Novaseq6010": "Illumina NovaSeq 6000",
            "Novaseq6006": "Illumina NovaSeq 6000",
            "Novaseq6005": "Illumina NovaSeq 6000",
            "Novaseq6004": "Illumina NovaSeq 6000",
            "Novaseq6003": "Illumina NovaSeq 6000",
            "Novaseq6002": "Illumina NovaSeq 6000",
            "Novaseq6001": "Illumina NovaSeq 6000",
            "Nextseq500-NS500488": "NextSeq 550",
        }
        instrument_model = instrument_model[
            metadata["ingest_metadata"]["metadata"]["acquisition_instrument_model"]
        ]

        assay_type = metadata["dataset_type"]

        # @icaoberg ignore field if missing. it is not known if field was moved
        if "preparation_protocol_doi" in metadata["ingest_metadata"]["metadata"]:
            protocols_io_doi = metadata["ingest_metadata"]["metadata"]["preparation_protocol_doi"]
        elif "protocols_io_doi" in metadata["ingest_metadata"]["metadata"]:
            protocols_io_doi = metadata["ingest_metadata"]["metadata"]["protocols_io_doi"]
        else:
            protocols_io_doi = "None"

        acquisition_instrument_vendor = metadata["ingest_metadata"]["metadata"][
            "acquisition_instrument_vendor"
        ]
        acquisition_instrument_model = metadata["ingest_metadata"]["metadata"][
            "acquisition_instrument_model"
        ]
        sequencing_reagent_kit_raw = metadata["ingest_metadata"]["metadata"][
            "sequencing_reagent_kit"
        sequencing_reagent_kit = sequencing_reagent_kit_raw.replace(';', '')
        ]

        # @icaoberg link is needed to map to a protocol description
        protocols_io = {
            "10.17504/protocols.io.86khzcw": "10X Genomics Single-Nucleus RNA-Sequencing for Transcriptomic Profiling of Adult Human Tissues V.3",
            "10.17504/protocols.io.bpgzmjx6": "Library Generation using Slide-seqV2 V.1",
            "10.17504/protocols.io.be5gjg3w": "SNARE-seq2 V.1",
            "10.17504/protocols.io.bfwajpae": "10X snRNAseq Nuclei Isolation and Library Preparation Protocol ",
            "10.17504/protocols.io.6t8herw": "Isolation of nuclei from frozen tissue for ATAC-seq and other epigenomic assays V.1",
            "10.17504/protocols.io.bf33jqqn": "10x Single Cell ATACseq Nuclei Isolation and Library Preparation Protocol",
            "10.17504/protocols.io.bvbmn2k6": "Chromium Multiome ATAC_GEX (10x)",
            "10.17504/protocols.io.bfsmjnc6": "Bulk WGS - Ultra II DNA Library Prep Kit for Illumina E7645/E7103",
            "10.17504/protocols.io.bfsnjnde": "BulkATAC - Isolation of nuclei from frozen tissue for ATAC-seq and other epigenomic assays",
            "10.17504/protocols.io.bftnjnme": "Bulk RNA - Protocol for use with NEBNext Poly(A) mRNA Magnetic Isolation Module (NEB #E7490) and NEBNext Ultra II Directional RNA Library Prep Kit for Illumina (E7760, E7765)",
            "10.17504/protocols.io.bukqnuvw": "Nuclei Isolation from Tissue for 10x Multiome",
            "10.17504/protocols.io.be79jhr6": "HuBMAP UF TMC - 10x Genomics scRNAseq Modality Overview",
            "10.17504/protocols.io.9yih7ue": "sci-RNA-seq3",
            "10.17504/protocols.io.be8mjhu6": "sci-ATAC-seq3",
        }

        try:
             if "preparation_protocol_doi" in metadata["ingest_metadata"]:
                protocols_io_title = protocols_io[
                    metadata["ingest_metadata"]["metadata"]["preparation_protocol_doi"]
            else:
                protocols_io_title = protocols_io[
                    metadata["ingest_metadata"]["metadata"]["protocols_io_doi"]
            ]
        except:
            protocols_io_title = None

        # deprecated design_description(s)
        design_description = f"The protocol and materials for the {assay_type} library construction process can be found in the following protocols.io protocol: dx.doi.org/{protocols_io_doi}. The library was sequenced on the {acquisition_instrument_vendor} {acquisition_instrument_model} system using the {sequencing_reagent_kit} kit."
        design_description = f"The {assay_type} library was sequenced on the {acquisition_instrument_vendor} {acquisition_instrument_model} system using the {sequencing_reagent_kit} kit."
        design_description = f"“A full description of the protocol and materials used in the {assay_type} library construction process can be found on protocols.io under the following protocol - Add Protocol Title Here.”"

        # current design_description
        if protocols_io_title is None:
            design_description = f"The {assay_type} library was sequenced on the {acquisition_instrument_vendor} {acquisition_instrument_model} system using the {sequencing_reagent_kit} kit. A full description of the protocol and materials used in the {assay_type} library construction process."
        else:
            design_description = f"The {assay_type} library was sequenced on the {acquisition_instrument_vendor} {acquisition_instrument_model} system using the {sequencing_reagent_kit} kit. A full description of the protocol and materials used in the {assay_type} library construction process can be found on protocols.io under the following protocol - {protocols_io_title}."

        reference_genome_assembly = None
        alignment_software = None

        dataset = hubmapinventory.get(hubmap_id, token=token)

        if dataset.empty:
            print(f"Dataset {hubmap_id} has no inventory")
            return df

        dataset = dataset.sort_values("filename")
        dataset = dataset[
            (dataset["filename"].str.contains("fq.gz"))
            | (dataset["extension"] == ".fastq.gz")
        ]
        dataset = dataset[["filename", "md5"]]

        datum = [
            dbgap_study_id,
            hubmap_id,
            library_id,
            title,
            library_strategy,
            library_source,
            library_selection,
            library_layout,
            platform,
            instrument_model,
            design_description,
            reference_genome_assembly,
            alignment_software,
        ]

        bash_script_file = f"{dbgap_study_id}/script.sh"
        if Path(bash_script_file).exists():
            Path(bash_script_file).unlink()

        __print_to_file(bash_script_file, "#!/bin/bash\n\n")
        commands = ""
        for index, row in dataset.iterrows():
            if prepend_sample_id:
                datum.extend(["fastq", f'{hubmap_id}-{row["filename"]}', row["md5"]])
                commands = f'{commands}cp -v "{dataset_directory}/{row["filename"]}" {hubmap_id}-{row["filename"]}\n'
            else:
                datum.extend(["fastq", row["filename"], row["md5"]])
                commands = f'{commands}cp -v "{dataset_directory}/{row["filename"]}" {row["filename"]}\n'

        __print_to_file(bash_script_file, commands)
        data.append(datum)

    print("Creating dataframe")
    df = pd.DataFrame(data)

    print("Writing dataframe to Excel spreadsheet")
    output_file = f"{dbgap_study_id}/spreadsheet.xlsx"
    write_excel(df, output_file)

    output_file = f"{dbgap_study_id}-{datetime.today().strftime('%Y%m%d')}.zip"
    print(f"Compressing folder {dbgap_study_id} to {output_file}")
    if Path(output_file).exists():
        Path(output_file).unlink()

    __compress_folder(dbgap_study_id, output_file)

    return df


##########################################################################
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
    donor = donor.reindex(columns=["SUBJECT_ID", "CONSENT", "SEX"])
    donor.to_csv(f"{directory}/2a_SubjectConsent_DS.txt", index=False, sep="\t")


##########################################################################
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
            datum["sample_id"], token=token, instance="prod"
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


##########################################################################
def __create_sample_mapping(df: pd.DataFrame, token: str, directory: str):
    sample_mapping = df[["donor_hubmap_id", "sample_id"]]
    sample_mapping = sample_mapping.rename(
        columns={"donor_hubmap_id": "SUBJECT_ID", "sample_id": "SAMPLE_ID"}
    )
    sample_mapping.to_csv(f"{directory}/3a_SSM_DS.txt", index=False, sep="\t")


##########################################################################
def __get_spreadhsheets(directory: str):
    filename = "2b_SubjectConsent_DD.xlsx"
    URL = "https://github.com/hubmapconsortium/py-hubmap-dbgap/blob/master/files/2b_SubjectConsent_DD.xlsx"
    temp_file = Path(f"{directory}/{filename}")
    response = requests.get(URL)
    temp_file.write_bytes(response.content)

    filename = "3b_SSM_DD.xlsx"
    URL = "https://github.com/hubmapconsortium/py-hubmap-dbgap/blob/master/files/3b_SSM_DD.xlsx"
    temp_file = Path(f"{directory}/{filename}")
    response = requests.get(URL)
    temp_file.write_bytes(response.content)

    filename = "6b_SampleAttributes_DD.xlsx"
    URL = "https://github.com/hubmapconsortium/py-hubmap-dbgap/blob/master/files/6b_SampleAttributes_DD.xlsx"
    temp_file = Path(f"{directory}/{filename}")
    response = requests.get(URL)
    temp_file.write_bytes(response.content)

    filename = "spreadsheet.xlsx"
    URL = "https://github.com/hubmapconsortium/py-hubmap-dbgap/raw/master/files/spreadsheet.xlsx"
    temp_file = Path(f"{directory}/{filename}")
    response = requests.get(URL)
    temp_file.write_bytes(response.content)


##########################################################################
def write_excel(
    df: pd.DataFrame, output_file: str
) -> None:  # Drops any completely empty rows/columns.
    df.dropna(how="all", inplace=True)
    # Adds non-file columns (indices 0-12, columns A-M).
    with pd.ExcelWriter(
        output_file, mode="a", if_sheet_exists="overlay", engine="openpyxl"
    ) as writer:
        df.to_excel(
            writer,
            sheet_name="Sequence_Data",
            startrow=1,
            startcol=0,
            index=False,
            header=False,
        )

    """
    Creates the repeated set of headers for each file, based on the number
    of columns minus 19 (representing spreadsheet columns A-M and
    the six file columns that exist already, N-S) divided by 3.
    It then appends that list to the end of the first row.
    """
    unique_files = int((len(df.columns) - 19) / 3)
    cols = unique_files * ["filetype", "filename", "MD5_checksum"]
    column_names = pd.DataFrame(cols)
    with pd.ExcelWriter(
        output_file, mode="a", if_sheet_exists="overlay", engine="openpyxl"
    ) as writer:
        column_names.T.to_excel(
            writer,
            sheet_name="Sequence_Data",
            startrow=0,
            startcol=19,
            index=False,
            header=False,
        )
