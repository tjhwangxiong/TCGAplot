% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gene_gender.R
\name{gene_gender}
\alias{gene_gender}
\title{gene_gender}
\usage{
gene_gender(cancer, gene, add = "jitter", palette = "jco", legend = "none")
}
\arguments{
\item{cancer}{cancer name likes "BLCA".}

\item{gene}{gene name likes "KLF7".}

\item{add}{character vector for adding another plot element. likes "none", "dotplot", "jitter".}

\item{palette}{the color palette to be used for coloring or filling by groups. Allowed values include scientific journal palettes from ggsci R package, e.g.: "npg", "aaas", "lancet", "jco", "ucscgb", "uchicago", "simpsons" and "rickandmorty".}

\item{legend}{legend position. Allowed values include "top","bottom","left","right" and "none".}
}
\description{
Create a box plot for a single gene with symbols indicating statistical significance grouped by gender in a specific type of cancer.
}
\examples{
gene_gender("BLCA","KLF7")
}
