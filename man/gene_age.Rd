% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gene_age.R
\name{gene_age}
\alias{gene_age}
\title{gene_age}
\usage{
gene_age(
  cancer,
  gene,
  age = 60,
  add = "jitter",
  palette = "jco",
  legend = "none"
)
}
\arguments{
\item{cancer}{cancer name likes "ACC".}

\item{gene}{gene name likes "KLF7".}

\item{age}{numeric number of age like 60.}

\item{add}{character vector for adding another plot element. likes "none", "dotplot", "jitter".}

\item{palette}{the color palette to be used for coloring or filling by groups. Allowed values include scientific journal palettes from ggsci R package, e.g.: "npg", "aaas", "lancet", "jco".}

\item{legend}{legend position. Allowed values include "top","bottom","left","right" and "none".}
}
\description{
Create a box plot for a single gene with symbols indicating statistical significance grouped by age in a specific type of cancer.
}
\examples{
gene_age("ACC","KLF7")
}
