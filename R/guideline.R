# Guideline for writing new `read_xxx` functions
#
# This document provides a guideline for writing new `read_xxx` functions in the `glyread` package.
# When writing a new `read_xxx` function, please follow the steps below:
# 1. Implement a reading function named `.read_xxx_df()`.
#    It should just read in the raw data file with column types specified.
# 2. Implement a tidying function named `.tidy_xxx()`.
#    It should take the output of `.read_xxx_df()` and perform the following steps:
#    1. Rename all neccessary columns into standard names
#       ("peptide", "peptide_site", "protein", "protein_site", "gene", "glycan composition", "glycan_structure").
#    2. Perform column transformation. Note that composition and structure parsing is not performed here.
#    3. Add neccesary columns.
#    4. If the raw data is in wide format, pivot it to long format.
#       The final `tidy_df` should have the columns mention above, along with "sample" and "value" columns.
# 3. In the API function, call `read_xxx_df()`, `.tidy_xxx()`, and `.read_template()` sequentially.
#    `.read_template()` takes care of all rest steps including aggregation and others.