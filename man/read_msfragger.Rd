% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/msfragger.R
\name{read_msfragger}
\alias{read_msfragger}
\title{Read MSFragger-Glyco result}
\usage{
read_msfragger(
  fp,
  sample_info = NULL,
  quant_method = c("label-free", "TMT"),
  glycan_type = c("N", "O"),
  sample_name_converter = NULL
)
}
\arguments{
\item{fp}{File path of the MSFragger result file.}

\item{sample_info}{File path of the sample information file (csv),
or a sample information data.frame/tibble.}

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
MSFragger-Glyco is a software for glycopeptide identification and quantification.
This function reads in the result file and returns a \code{\link[glyexp:experiment]{glyexp::experiment()}} object.

Because MSFragger performs quantification on the PSM level for each sample,
and outputs the result to separate files,
this function cannot read and merge multiple result files into a single \code{\link[glyexp:experiment]{glyexp::experiment()}} object.
See details to understand why.

As a workaround, use this code snippet to read the result into a list of \verb{glyexp::experiment()`` objects, aggregate the experiments to the "glycoform" level using }glyclean::aggregate()\verb{, and then use }glyexp::merge()` (need glyexp >= 0.7.0) to merge them into one:

\if{html}{\out{<div class="sourceCode r">}}\preformatted{library(tidyverse)
library(glyread)
library(glyclean)
library(glyexp)

# Read the result into a list of experiments
files <- c("psm1.tsv", "psm2.tsv", "psm3.tsv")
exps <- map(files, read_msfragger)  # ignore sample_info for now

# Aggregate the experiments to the "glycoform" level
exps_agg <- map(exps, aggregate, to_level = "gf")  # or other levels

# Merge the experiments into one
exp_merged <- reduce(exps_agg, merge)

# Add sample information
# Say you already have a sample information tibble (`sample_df`) with a "sample" column
exp_merged$sample_info <- exp_merged$sample_info |>
  left_join(sample_df, by = "sample")
}\if{html}{\out{</div>}}

Note that we override the sample information with the new one,
which is not a recommended practice.
However, the code here is harmless,
and you can safely use it until we find a better solution.
}
\section{Why not read multiple files?}{
This is the first solution we thought of.
However, we found it is not compatible with the current design of \code{glycoverse} settings.
MSFragger-Glyco results are on the PSM level,
each raw sample resulting in a separate file.
To merge them into one expression matrix, with rows as PSMs and columns as samples,
we need to identify the PSM that appears in different samples.
This is not trivial, as an ion (with the same charge and modification state)
can generate multiple PSMs with different retention time and m/z.

One possible solution is to aggregate each PSM result into a glycoform level.
That is, for a glycoform (unique combination of glycan, protein, and glycosite),
we sum up the quantifications of all PSMs that belong to this glycoform.
Then the data can be merged into one expression matrix,
with rows as glycoforms and columns as samples.
This is exactly what the about code snippet does.

So why not just incorporate this approach into \code{glyread}?
Well, the job of \code{glyread} should be to read search results as they are.
It should not be responsible for any data preprocessing.
Aggregating data into "glycoform" level is a preprocessing step.
In fact, it is not even the first step
(aggregation after normalization and imputation is more preferred).
So we put this function in the \code{glyclean} package.

Maybe in the future, MSFragger-Glyco has a merged output like pGlycoQuant.
Then everyone is happy.
But for now, we leave this to the user to do.
Hint: you can wrap the snippet above into a function yourself.
}

\seealso{
\code{\link[glyexp:experiment]{glyexp::experiment()}}, \code{\link[glyrepr:glycan_composition]{glyrepr::glycan_composition()}}
}
