% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gene_coexp_heatmap.R
\name{gene_coexp_heatmap}
\alias{gene_coexp_heatmap}
\title{gene_coexp_heatmap}
\usage{
gene_coexp_heatmap(cancer, gene, top_n = 20, method = "pearson")
}
\arguments{
\item{cancer}{cancer name likes "STAD".}

\item{gene}{gene name likes "KLF7".}

\item{top_n}{the number of co-expressed genes.}

\item{method}{method="pearson" is the default value. The alternatives to be passed to correlation were "spearman" and "kendall".}
}
\description{
Heatmap and Go enrichment of the positive and negative co-expressed genes of a single gene in a specific type of cancer.
}
\examples{
gene_coexp_heatmap("STAD","KLF7")
}
