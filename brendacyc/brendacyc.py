"""Main module."""
from pathlib import Path
from typing import List

import pandas as pd


class BrendaDB:
    def __init__(self, filepath: str) -> None:
        self.fp = Path(filepath).expanduser().resolve()
        if not self.fp.is_file():
            raise FileNotFoundError(self.fp)

        self.brenda_fields = {
            "ACTIVATING_COMPOUND",
            "APPLICATION",
            "CLONED",
            "COFACTOR",
            "CRYSTALLIZATION",
            "ENGINEERING",
            "EXPRESSION",
            "GENERAL_INFORMATION",
            "GENERAL_STABILITY",
            "IC50_VALUE",
            "INHIBITORS",
            "KCAT_KM_VALUE",
            "KI_VALUE",
            "KM_VALUE",
            "LOCALIZATION",
            "METALS_IONS",
            "MOLECULAR_WEIGHT",
            "NATURAL_SUBSTRATE_PRODUCT",
            "ORGANIC_SOLVENT_STABILITY",
            "OXIDATION_STABILITY",
            "PH_OPTIMUM",
            "PH_RANGE",
            "PH_STABILITY",
            "PI_VALUE",
            "POSTTRANSLATIONAL_MODIFICATION",
            "PROTEIN",
            "PURIFICATION",
            "REACTION",
            "REACTION_TYPE",
            "RECOMMENDED_NAME",
            "REFERENCE",
            "RENATURED",
            "SOURCE_TISSUE",
            "SPECIFIC_ACTIVITY",
            "STORAGE_STABILITY",
            "SUBSTRATE_PRODUCT",
            "SUBUNITS",
            "SYNONYMS",
            "SYSTEMATIC_NAME",
            "TEMPERATURE_OPTIMUM",
            "TEMPERATURE_RANGE",
            "TEMPERATURE_STABILITY",
            "TRANSFERRED_DELETED",
            "TURNOVER_NUMBER",
        }

    def read_brenda(self, clean: bool = True) -> pd.DataFrame:
        """Reads raw BRENDA text file and load all non-empty lines.

        Empty lines and comment lines (starting with `*`) are skipped.
        The text file should be downloaded from:
        https://www.brenda-enzymes.org/download_brenda_without_registration.php

        Args:
            filepath (str): A string indicating the path to the text file.

        Returns:
            List[str]: Each element is a line in the file.
        """
        with open(self.fp, "r") as f:
            brenda_txt = []
            for line in f:
                line = line.rstrip()
                if line and line[0] != "*":
                    brenda_txt.append(line)

        df = self._txt2df(brenda_txt)
        if clean:
            df = self._clean_ec_number(df)
        return df

    def _txt2df(self, brenda_txt: List[str]) -> pd.DataFrame:
        """Converts list of strings to a data frame.

        Args:
            brenda_txt (List[str]): Lines extracted from a BRENDA text file.

        Returns:
            pd.DataFrame: For each EC entry, the annotations are split into three columns:
            - ID: EC number, e.g. 1.1.1.1
            - field: the content of the information, one of `self.brenda_fields`
            - description: everything else
        """
        colID, colField, colDescription = [], [], []

        current_ID = brenda_txt[0][3:]  # ID\tx.x.x.x, remove ID\t
        current_field = brenda_txt[1]  # one of self.brenda_fields
        ec_info = ""

        i = 2  # Skip first two lines because they're already read
        nline = len(brenda_txt) - 1  # Not reading last line
        while i < nline:
            if brenda_txt[i] == "///":
                # /// indicates the end of an EC-number specific part, so we should
                # insert the previous entry
                colID.append(current_ID)
                colField.append(current_field)
                colDescription.append(ec_info)

                # Update ID and first field, then clear ec_info to prepare for the
                # next chunk
                i += 1
                current_ID = brenda_txt[i][3:]  # drop the leading ID\t
                i += 1
                current_field = brenda_txt[i]
                ec_info = ""
            else:
                if brenda_txt[i] in self.brenda_fields:
                    # Insert previous entry, update field and clear ec_info
                    # when we hit the next field
                    colID.append(current_ID)
                    colField.append(current_field)
                    colDescription.append(ec_info)
                    current_field = brenda_txt[i]
                    ec_info = ""
                else:
                    # Append line to the current field as a continuation
                    # of the previous lines
                    ec_info = ec_info + brenda_txt[i] + "\n"

            i += 1

        # Insert last entry
        colID.append(current_ID)
        colField.append(current_field)
        colDescription.append(ec_info)

        return pd.DataFrame(
            {"ID": colID, "field": colField, "description": colDescription}
        )

    def _clean_ec_number(self, df: pd.DataFrame) -> pd.DataFrame:
        """Handles deleted and transferred EC numbers.

        In the ID column, some EC numbers have comments wrapped in parentheses
        about deleted or transferred EC numbers. These entries are often duplicated
        with varying `field`s and empty `description`s. Thus we drop these rows
        and append a new row with the EC number as the `ID`, "TRANSFERRED_DELETED"
        as the `field`, and the comment in the `description` column.

        Empty comments are removed from the ID column, leaving just the EC numbers.
        Comments for deleted EC numbers are usually about why they were deleted.
        For transferred EC numbers, the comments point to the new EC numbers.

        Args:
            df (pd.DataFrame): Output of `self.read_brenda`.

        Returns:
            pd.DataFrame: Same data with cleaned EC numbers.
        """
        # Remove empty comments
        df["ID"] = df["ID"].str.replace(" ()", "", regex=False)

        # Separate out deleted and transferred entries
        df_standard = df[~df["ID"].str.contains("(", regex=False)]
        df_nonstd = (
            df.query("ID.str.contains('(', regex=False)")
            .drop_duplicates(subset=["ID"])
            .assign(field="TRANSFERRED_DELETED")
            .assign(description=lambda x: x["ID"].str.extract(r"\((.*)\)$"))
            .assign(ID=lambda x: x["ID"].str.replace(r"\s?\(.*$", "", regex=True))
        )

        # Append these entries at the end of the data frame
        return pd.concat([df_standard, df_nonstd], axis=0)


if __name__ == "__main__":
    db = BrendaDB("~/Downloads/brenda_download.txt")
    df = db.read_brenda()

