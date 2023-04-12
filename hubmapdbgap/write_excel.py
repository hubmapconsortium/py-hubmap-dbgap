import os
import pandas as pd

# Note: placeholder function name
def write_excel(df: pd.DataFrame, output_file: str | os.PathLike) -> None:
    # Drops any completely empty rows/columns.
    df.dropna(how="all", inplace=True)
    # Adds phs_accession and sample_IDm columns (indices 0 & 1).
    with pd.ExcelWriter(output_file, mode="a", if_sheet_exists="overlay") as writer:
        df.to_excel(
            writer,
            sheet_name="Sequence_Data",
            columns=df.columns[:2],
            startrow=1,
            startcol=0,
            index=False,
            header=False,
        )
    # Adds file data, beginning at index 2.
    with pd.ExcelWriter(output_file, mode="a", if_sheet_exists="overlay") as writer:
        df.to_excel(
            writer,
            sheet_name="Sequence_Data",
            columns=df.columns[2:],
            startrow=1,
            startcol=13,
            index=False,
            header=False,
        )
    """
    Creates the repeated set of headers for each file, based on the number
    of columns minus 8 (representing phs_accession, sample_IDm, and
    the six file columns that exist already) divided by 3.
    It then appends that list to the end of the first row.
    """
    cols = []
    unique_files = int((len(df.columns) - 8) / 3)
    for _ in range(unique_files):
        cols.extend(["filetype", "filename", "MD5_checksum"])
    column_names = pd.DataFrame(cols)
    with pd.ExcelWriter(output_file, mode="a", if_sheet_exists="overlay") as writer:
        column_names.T.to_excel(
            writer,
            sheet_name="Sequence_Data",
            startrow=0,
            startcol=19,
            index=False,
            header=False,
        )
