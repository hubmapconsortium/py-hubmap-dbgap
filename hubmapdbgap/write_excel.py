import os
import pandas as pd


def write_excel(df: pd.DataFrame, output_file: str | os.PathLike) -> None:
    # Drops any completely empty rows/columns.
    df.dropna(how="all", inplace=True)
    # Adds non-file columns (indices 0-12, columns A-M).
    with pd.ExcelWriter(output_file, mode="a", if_sheet_exists="overlay") as writer:
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
    with pd.ExcelWriter(output_file, mode="a", if_sheet_exists="overlay") as writer:
        column_names.T.to_excel(
            writer,
            sheet_name="Sequence_Data",
            startrow=0,
            startcol=19,
            index=False,
            header=False,
        )
