# GWAS QC & LiftOver Pipeline

A Python pipeline designed for automated Quality Control (QC), coordinate conversion (LiftOver), and formatting of GWAS summary statistics. This tool processes raw summary statistics (hg38), cleans them, converts them to hg19, filters specific regions, and outputs data ready for downstream analyses (e.g., FUMA, SMR, Coloc).

## 📋 Input File Format

The pipeline expects input files (e.g., `.txt` or `.tsv`) to be **space or tab-delimited**. 
The file **must** contain the following columns with these specific headers:

| Column Header | Description | Data Type | Example |
| :--- | :--- | :--- | :--- |
| **SNP** | Variant Identifier (rsID) | String | `rs12345` |
| **chr** | Chromosome number | Integer/String | `1` or `chr1` |
| **pos** | Base pair position (hg38) | Integer | `10500` |
| **effect_allele** | The effect allele (A1) | String | `A` |
| **other_allele** | The non-effect/other allele (A2) | String | `G` |
| **eaf** | Effect Allele Frequency | Float (0-1) | `0.15` |
| **beta** | Effect size (Beta coefficient) | Float | `0.02` |
| **se** | Standard Error | Float | `0.005` |
| **pval** | P-value | Float | `1.5e-08` |
| **samplesize** | Sample size (N) | Integer | `100000` |

> **Note:** Rows with missing values (`NaN`) in any of these critical columns will be automatically removed during the Deep Cleaning step.

---

## 🚀 Key Features & QC Pipeline

This pipeline performs the following steps sequentially:

* **Deep Cleaning**
    * Automatically validates data integrity.
    * Removes rows with missing values (`NaN`) or formatting errors in critical columns (e.g., `chr`, `pos`, `beta`, `pval`, `eaf`).
    * Standardizes chromosome notation.

* **Coordinate Conversion (LiftOver)**
    * Converts genomic coordinates from **hg38 to hg19** using `pyliftover`.
    * Automatically drops variants that fail conversion (e.g., variants located in gap regions).

* **MHC Region Filtering**
    * Removes potentially confounded variants located in the Major Histocompatibility Complex (MHC) region.
    * **Region defined (hg19):** `chr6:28,477,797–33,448,354`.

* **MAF Filtering**
    * Calculates Minor Allele Frequency (MAF) using the formula: `min(eaf, 1-eaf)`.
    * Removes rare variants with **MAF < 0.01**.

* **FUMA Formatting**
    * Generates compressed files (`.txt.gz`) strictly formatted for **FUMA** upload.
    * **Mapped Columns:** `SNP`, `A1`, `A2`, `BETA`, `SE`, `P`, `N`, `CHR`, `BP`.

* **Dual Version Output**
    * Saves the final cleaned datasets in both coordinate systems:
        1.  **hg19**: Cleaned and filtered data (ready for Magma/FUMA).
        2.  **hg38**: The cleaned hg19 data converted back to hg38 (consistent filtering).

## 📂 Output Structure

After running the script, the following directories will be created:

* `1.hg19/`: Cleaned summary statistics in hg19 coordinates.
* `2.hg38/`: Cleaned summary statistics in hg38 coordinates.
* `3.FUMA/`: `.txt.gz` files formatted for FUMA platform.

## 🛠 Dependencies

* Python 3.x
* pandas
* pyliftover
* tqdm
* numpy
