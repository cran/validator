\name{intVal}
\alias{intVal}

\title{
Internal Validation Indices
}
\description{
This function is calculating the values of certain internal validation indices.
}
\usage{
intVal(y, x, index = "all")
}

\arguments{
  \item{y}{
Object of class \code{kcca} returned by clustering methods from the package \code{flexclust}
}
  \item{x}{
Data matrix, which contains the observations of clustering or data matrix of a test data set
}
  \item{index}{
The internal validation indices, which are calculated: "\code{calinski}", "\code{db}", "\code{hartigan}", "\code{ratkowsky}", "\code{scott}", "\code{marriot}", "\code{ball}", "\code{trcovw}", "\code{tracew}", "\code{friedman}", "\code{rubin}", "\code{xuindex}", "\code{dunn}", "\code{connectivity}", "\code{silhouette}", and "\code{all}".
}
}
\details{
The internal validation indices based all on the references below. The indices are defined as:

\emph{calinski}: \eqn{\frac{SSB/(k-1)}{SSW/(n-k)}}, where \eqn{SSB} is the sum of squares between, 
\eqn{SSW} is the sum of squares between, \eqn{n} is the number of data points, 
and \eqn{k} is the number of clusters.

\emph{db}: \eqn{\frac{1}{k} \sum_{i=1}^k R_i}, where \eqn{R_i} as the maximum value of \eqn{R_{ij}=\frac{s_i + s_j}{d_{ij}}} with \eqn{s_i} is the similarity in the clusters and \eqn{d_{ij}} is the dissimilarity between the two clusters.

\emph{hartigan}: \eqn{\log \frac{SSB}{SSW}}

\emph{ratkowsky}: \eqn{\tilde{c}/ \sqrt{k}}, where \eqn{\tilde{c} = mean
  \sqrt{varSSB / varSST}} and the abbreviation emph{var} stands for each variable and \eqn{SST} for the sum of squares total.

\emph{ball}: \eqn{\frac{SSW}{k}}

\emph{xuindex}: \eqn{d\log(\sqrt{SSW/(dn^2)}) + \log(k)}, with \eqn{d} as the dimension of the data points.

\emph{scott}: \eqn{n \log \frac{\det(T)}{\det(W)}}, with \eqn{T} is the scatter distance matrix and \eqn{W} is the pooled within groups scatter matrix.

\emph{marriot}: \eqn{k^2 \det(W)}

\emph{trcovw}: \eqn{trace(cov W)}

\emph{tracew}: \eqn{trace(W)}

\emph{friedman}: \eqn{trace(W^{-1}B)}, where \eqn{B} is the between groups scatter matrix.

\emph{rubin}: \eqn{\det(T) / \det(W)}

The R code of these functions above are taken from the package \code{cclust} with small changes, e.g. by computing the \eqn{SSW} regarding their cluster centers.

\emph{dunn}: \eqn{D(\mathcal{C}) = \frac{\min\limits_{C_k, C_l \in \mathcal{C}, C_k \neq
    C_l} \Big( \min\limits_{x_i \in C_k, x_j \in C_l} dist(x_i, x_j) \Big)}
{\max\limits_{C_m \in \mathcal{C}} diam (C_m)}}, where \eqn{diam(C_m)} is maximum distance between each data item in the
cluster \eqn{C_m}, and \eqn{dist(x_i, x_j)} is the distance between the pairs of the
data points.

\emph{silhouette}: \eqn{S (x_i) = \frac{b_i - a_i}{max(a_i, b_i)}}, where \eqn{a_i} is the average distance between the data points to all the other observations in the same cluster and \eqn{b_i} is the average distance between data points to all other points from the closest neighbouring cluster. Then, the average of all Silhouette Widths is computed.

\emph{connectivity}: \eqn{Conn( \mathcal{C} ) = \sum_{i=1}^{N} \sum_{j=1}^{L} x_{i, nn_{i(j)}}}, where \eqn{nn_{i(j)}} be the j-th neighbour of the data point \eqn{x_i }; so \eqn{x_{i, nn_{i(j)}}} is zero if \eqn{x_i} and \eqn{nn_{i(j)}} are in
the same cluster, otherwise the value is computed by \eqn{\frac{1}{j}} and \eqn{L} determines the number
of neighbours that contribute to the connectivity measure.


}
\value{
This function returns a vector with the internal validation indices.
}
\references{
Glenn W. Milligan and Martha C. Cooper, \emph{An examination of procedures for de- termining the number of clusters in a dataset}.

Julia Handl, Joshua Knowles, and Douglas B. Kell, \emph{Computational cluster validation in post-genomic data analysis},
\url{http://dbkgroup.org/handl/clustervalidation/}.

Andreas Weingessel, Evgenia Dimitriadou, and Sara Dolnicar, \emph{An examination of indexes for determining the number of clusters in binary data sets}.

Guy Brock, Vasyl Pihur, Susmita Datta, and Somnath Datta, \emph{clValid: An R package for cluster validation},
\url{http://www.jstatsoft.org/v25/i04}.
}
\author{
Marcus Scherl
}

\examples{
# require(mlbench)
# require(flexclust)

x <- mlbench.2dnormals(500, 3)
cl <- kcca(x$x, 3)
intVal(cl, x$x)

x <- mlbench.2dnormals(500, 3)
cl <- kmeans(x$x, 3)
cl <- as.kcca(cl, x$x)
intVal(cl, x$x)

}
\keyword{intVal}