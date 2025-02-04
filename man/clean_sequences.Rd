% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clean_sequences.R
\name{clean_sequences}
\alias{clean_sequences}
\title{Remove non-standard characters from sequences}
\usage{
clean_sequences(sequences, type, alphabet = NULL, replacement_character = NULL)
}
\arguments{
\item{sequences}{A character vector of sequences}

\item{type}{"aa" or "nt"}

\item{alphabet}{(optional) Provide a custom alphabet (defaults for "aa" and "nt" are the alphabet or "-")}

\item{replacement_character}{(optional) Provide a custom replacement character (defaults for "aa" and "nt" are "X" and "N"; you may want "-")}
}
\value{
a character vector of cleaned sequences
}
\description{
Remove non-standard characters from sequences
}
