% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pglyco3-pglycoquant.R
\name{read_pglyco3_pglycoquant}
\alias{read_pglyco3_pglycoquant}
\title{Read pGlyco3-pGlycoQuant result}
\usage{
read_pglyco3_pglycoquant(
  fp,
  sample_info = NULL,
  name = NULL,
  quant_method = c("label-free", "TMT"),
  glycan_type = c("N", "O"),
  sample_name_converter = NULL
)
}
\arguments{
\item{fp}{File path of the pGlyco3 result file.}

\item{sample_info}{File path of the sample information file,
or a sample information tibble.}

\item{name}{Name of the experiment. If not provided, a default name with
current time will be used.}

\item{quant_method}{Quantification method. Either "label-free" or "TMT".}

\item{glycan_type}{Glycan type. Either "N" or "O". Default is "N".}

\item{sample_name_converter}{A function to convert sample names from file paths.
The function should take a character vector of old sample names
and return new sample names.
Note that sample names in \code{sample_info} should match the new names.
If NULL, original names are kept.}
}
\value{
An \code{\link[glyexp:experiment]{glyexp::experiment()}} object.
}
\description{
If you used pGlyco3 for intact glycopeptide identification,
and used pGlycoQuant for quantification, this is the function for you.
It reads in a pGlycoQuant result file and returns a \code{\link[glyexp:experiment]{glyexp::experiment()}} object.
Currently only label-free quantification is supported.
}
\section{Input}{
You should use the "Quant.spectra.list" file in the pGlycoQuant result folder.
Files from pGlyco3 result folder are not needed.
For instructions on how to use pGlyco3 and pGlycoQuant, please refer to
the manual: \href{https://github.com/Power-Quant/pGlycoQuant/blob/main/Manual\%20for\%20pGlycoQuant_v202211.pdf}{pGlycoQuant}.

The sample information file should be a \code{csv} file with the first column
named \code{sample}, and the rest of the columns being sample information.
The \code{sample} column must match the \code{RawName} column in the pGlyco3 result file,
although the order can be different.

You can put any useful information in the sample information file.
Recommended columns are:
\itemize{
\item \code{group}: grouping or conditions, e.g. "control" or "tumor",
required for most downstream analyses
\item \code{batch}: batch information, required for batch effect correction
\item \code{bio_replicate}:  Indicates the identity of biologically distinct samples.
While the \code{sample} column is typically used as a unique identifier for each run,
\code{bio_replicate} specifies whether the samples originate from the same biological source.
In common experimental designs,
multiple technical replicates are often derived from the same biological sample,
making the bio_replicate column essential for distinguishing
biological replicates from technical replicates.
}
}

\section{Output}{
This function returns a \code{\link[glyexp:experiment]{glyexp::experiment()}} object.

The following columns could be found in the variable information tibble:
\itemize{
\item \code{charge}: integer, charge state
\item \code{peptide}: character, peptide sequence
\item \code{modifications}: character, modifications other than glycosylation,
separated by semicolon, e.g. \verb{5,Carbamidomethyl[C];10,Carbamidomethyl[C]}
\item \code{glycan_composition}: character, glycan composition, e.g. "H5N4F1A1"
\item \code{glycan_structure}: character, pGlyco-style structure strings, renamed from
\code{PlausibleStruct} in the original file
\item \code{peptide_site}: integer, site of glycosylation on peptide
\item \code{proteins}: character, protein names, separated by semicolon
\item \code{genes}: character, gene names, separated by semicolon
\item \code{protein_sites}: character, site of glycosylation on protein,
separated by semicolon
}
}

