% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methy_kmplot.R
\name{methy_kmplot}
\alias{methy_kmplot}
\title{methy_kmplot}
\usage{
methy_kmplot(cancer, gene, palette = "jco")
}
\arguments{
\item{cancer}{cancer name likes "COAD".}

\item{gene}{gene name likes "KLF7".}

\item{palette}{the color palette to be used for coloring or filling by groups. Allowed values include scientific journal palettes from ggsci R package, e.g.: "npg", "aaas", "lancet", "jco".}
}
\description{
Describes the K_M survival plot based on the promoter methylation of a single gene in a specific type of cancer. A pdf file named methylation_kmplot will be generated in the working directory.
}
\examples{
methy_kmplot("COAD","KLF7")
}
