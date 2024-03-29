% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SPARKS_analysis.R
\name{calculate_fgsea_score}
\alias{calculate_fgsea_score}
\title{Calculate enrichment score using FGSEA}
\usage{
calculate_fgsea_score(gsea_library, study_rank)
}
\arguments{
\item{gsea_library}{Signature event list (positive and negative) from `extract_GSEA_significant_events` function}

\item{study_rank}{AS events sorted by delta PSI}
}
\value{
gsea_result_df Dataframe with enrichment score and pval for positive and negative event sets
}
\description{
Calculate enrichment score using FGSEA
}
