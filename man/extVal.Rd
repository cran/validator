\name{extVal}
\alias{extVal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
External Validation Indices
}
\description{
This function is calculating the values of certain external validation indices. The values compare the predicted cluster assignment with their true cluster labelling existing in the data.  
}
\usage{
extVal(x, y, index = "all")
}

\arguments{
  \item{x}{
Vector of the true cluster assignment
}
  \item{y}{
Vector of the cluster assignment, which has to be compare with the true one
}
  \item{index}{
The external indices, which are calculated: "\code{Hamann}", "\code{Czekanowski}", "\code{Kulczynski}", "\code{McConnaughey}", "\code{Peirce}", "\code{Wallace1}", "\code{Wallace2}", "\code{Gamma}", "\code{Sokal1}", "\code{Fager}", "\code{Sokal2}", "\code{Sokal3}", "\code{Gower}", "\code{Roger}", "\code{Kruskal}", "\code{Pearson}", "\code{Rand}", "\code{Jaccard}", "\code{Folkes}", "\code{Russel}", and "\code{all}".
}
}
\details{
The data points are counted in a pairwise co-assignment. Given two partitions named C1 and C2, the quantities a, b, c and d are computed for the pairs of the data points xi and xj and their cluster assignments. The numbers a and d are counted as the agreements between the two cluster partitions C1 and C2, whereas b and c are the disagreements of these two partitions:
\tabular{ll}{
		a \tab - number of pairs in the same cluster of 
		both partitions, \cr
		b \tab - number of pairs, which are in the same cluster of partition C1 but 
		not of partition C2, \cr
		c \tab - number of pairs, in which both data points are in different clusters in the partitions C1 and C2, \cr
		d \tab - number of pairs in different clusters of both partitions.
	}
Then, the external indices are computed with these values:
\tabular{ll}{
Rand \tab \eqn{\frac{a+d}{a+b+c+d}} \cr
Hamann \tab \eqn{\frac{(a+d)-(b+c)}{a+b+c+d}} \cr
Czekanowski \tab \eqn{\frac{2a}{2a+b+c}} \cr
Kulczynski \tab \eqn{\frac{1}{2} ( \frac{a}{a+b} + \frac{a}{a+c} ) } \cr
McConnaughey \tab \eqn{\frac{a^2-bc}{(a+b)(a+c)}} \cr
Peirce \tab \eqn{\frac{ad-bc}{(a+c)-(b+d)}} \cr
Fowlkes \tab \eqn{\frac{a}{\sqrt{(a+b)(a+c)}}} \cr
Wallace1 \tab \eqn{\frac{a}{a+b}}  \cr
Wallace2 \tab \eqn{\frac{a}{a+c}}  \cr
Gamma \tab \eqn{\frac{ad-bc}{\sqrt{(a+b)(a+c)(c+d)(b+d)}}} \cr
Sokal1 \tab \eqn{\frac{1}{4} ( \frac{a}{a+b} + \frac{a}{a+c} + \frac{d}{d+b} + \frac{d}{d+c})} \cr
Russel \tab \eqn{\frac{a}{a+b+c+d}} \cr
Fager \tab \eqn{\frac{a}{\sqrt{(a+b)(a+c)}} - \frac{1}{2 \sqrt{(a+b)}}} \cr
Pearson \tab \eqn{\frac{ad-bc}{(a+b)(a+c)(c+d)(b+d)}} \cr
Jaccard \tab \eqn{\frac{a}{a+b+c}} \cr
Sokal2 \tab \eqn{\frac{a}{a+2(b+c)}} \cr
Sokal3 \tab \eqn{\frac{ad}{\sqrt{(a+b)(a+c)(c+d)(b+d)}}} \cr
Gower \tab \eqn{\frac{a+d}{a+\frac{1}{2}(b+c)+d}} \cr
Rogers \tab \eqn{\frac{a+d}{a+2(b+c)+d}} \cr
Kruskal \tab \eqn{\frac{ad-bc}{ad+bc}} 

}
}
\value{
This function returns a vector with the external validation indices.
}
\references{
Ahmed N Albatineh, Magdalena Niewiadomska-Bugaj, and Daniel Mihalko \emph{On similarity indices and correction for chance agreement}.
}
\author{
Marcus Scherl
}

\examples{
#require(mlbench)
#require(flexclust)

x <- mlbench.2dnormals(500, 3)
cl <- kcca(x$x, 3)
pred <- predict(cl, x$x)
extVal(pred, x$class)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{extVal}
