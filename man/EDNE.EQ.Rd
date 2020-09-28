\name{EDNE.EQ}
\alias{EDNE.EQ}
\concept{multivariate statistics}
\concept{equivalence}
\concept{two-sample design}
\concept{Euclidean distance}

\title{
The EDNE-test for equivalence
}
\description{
The function \code{EDNE.EQ()} implements the EDNE-test for equivalence according to Hoffelder et al. (2015). It is a multivariate two-sample equivalence procedure. Distance measure of the test is the Euclidean distance.  
}
\usage{
EDNE.EQ(X, Y, eq_margin, alpha = 0.05, print.results = TRUE)
}
\arguments{
  \item{X}{
numeric data matrix of the first sample (REF). The rows of \code{X} contain the individual observations of the REF sample, the columns contain the variables/components of the multivariate sample. 
}
  \item{Y}{
numeric data matrix of the second sample (TEST). The rows of \code{Y} contain the individual observations of the TEST sample, the columns contain the variables/components of the multivariate sample. 
}
  \item{eq_margin}{
numeric (>0). The equivalence margin of the test. 
}
  \item{alpha}{
numeric (0<\code{alpha}<1). The significance level of the EDNE-test for equivalence. Usually set to 0.05 which is the default.
}
  \item{print.results}{
logical; if TRUE (default) summary statistics and test results are printed in the output. If NO no output is created 
}
}
\details{
This function implements the EDNE-test for equivalence. Distance measure of the test is the Euclidean distance. The test is an asymptotically valid test for the family of distributions fulfilling the assumptions of the multivariate central limit theorem (for further details see Hoffelder et al.,2015).   
}
\value{
a data frame; three columns containing the results of the test
  \item{p.value}{numeric; the p-value of the equivalence test according to Hoffelder et al. (2015)}
  \item{testresult.num}{numeric; 0 (null hypothesis of nonequivalence not rejected) or 1 (null hypothesis of nonequivalence rejected, decision in favor of equivalence)}
  \item{testresult.text}{character; test result of the test in text mode}
}
\references{
Hoffelder, T., Goessl, R., Wellek, S. (2015). Multivariate Equivalence Tests for Use in Pharmaceutical Development. \emph{Journal of Biopharmaceutical Statistics}, 25:3, 417-437.
URL: \url{http://dx.doi.org/10.1080/10543406.2014.920344}
}
\author{
Thomas Hoffelder <thomas.hoffelder at boehringer-ingelheim.com>
}

\examples{
# A recalculation of the three-dimensional EDNE example evaluation 
# in Hoffelder et al. (2015) can be done with the following code:

data(ex_data_JoBS)
REF_JoBS <- cbind(ex_data_JoBS[ which(ex_data_JoBS$Group=='REF'), ]
            [c("Diss_15_min","Diss_20_min","Diss_25_min")])
TEST_JoBS <- cbind(ex_data_JoBS[ which(ex_data_JoBS$Group=='TEST'), ]
            [c("Diss_15_min","Diss_20_min","Diss_25_min")])
equivalence_margin_JoBS <- 297
test_EDNE_JoBS <- EDNE.EQ(X=REF_JoBS
                          , Y=TEST_JoBS
                          , eq_margin=equivalence_margin_JoBS
                          , print.results = TRUE)
}

\keyword{multivariate}
\keyword{robust}