% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/import_data.R
\name{import_MATS}
\alias{import_MATS}
\title{Main import function to import MATS information into object}
\usage{
import_MATS(object, MATS_file, event_type)
}
\arguments{
\item{object}{the object}

\item{MATS_file}{MATS file from processing pipeline}

\item{event_type}{event_type}
}
\value{
object with the list of dataframe added
}
\description{
This function reads in the file and add it to the object
MATS result is included to perform differential splicing analysis
}
