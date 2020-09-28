\name{EDNE.EQ-package}
\alias{EDNE.EQ-package}
\concept{multivariate statistics}
\concept{equivalence}
\concept{two-sample design}
\concept{Euclidean distance}
\concept{dissolution profiles}

\docType{package}
\title{
\packageTitle{EDNE.EQ}
}
\description{
\packageDescription{EDNE.EQ}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{EDNE.EQ}
\packageIndices{EDNE.EQ}
}
\author{
\packageAuthor{EDNE.EQ}

Maintainer: \packageMaintainer{EDNE.EQ}
}
\references{
EMA (2010). Guidance on the Investigation of Bioequivalence.
European Medicines Agency, CHMP, London. Doc. Ref.: CPMP/EWP/QWP/1401/98 Rev. 1/ Corr **. 
URL: \url{https://www.ema.europa.eu/en/documents/scientific-guideline/guideline-investigation-bioequivalence-rev1_en.pdf}

FDA (1997). Guidance for Industry: Dissolution Testing of Immediate Release Solid Oral Dosage Forms. 
Food and Drug Administration FDA, CDER, Rockville.
URL: \url{https://www.fda.gov/media/70936/download}

Hoffelder, T., Goessl, R., Wellek, S. (2015). Multivariate Equivalence Tests for Use in Pharmaceutical Development. \emph{Journal of Biopharmaceutical Statistics}, 25:3, 417-437.
URL: \url{http://dx.doi.org/10.1080/10543406.2014.920344}

Suarez-Sharp, S., Abend, A., Hoffelder, T., Leblond, D., Delvadia, P., Kovacs, E., Diaz, D.A. (2020). In Vitro Dissolution Profiles Similarity Assessment in Support of Drug Product Quality: What, How, When - Workshop Summary Report. \emph{The AAPS Journal}, 22:74.
URL: \url{http://dx.doi.org/10.1208/s12248-020-00458-9}
}

\keyword{ package }
\keyword{multivariate}
\keyword{robust}

\examples{
# A recalculation of the three-dimensional EDNE example evaluation 
# in Hoffelder et al. (2015) can be done with the following code:

data(ex_data_JoBS)
REF_JoBS <- cbind(ex_data_JoBS[ which(ex_data_JoBS$Group=='REF'), ]
                  [c("Diss_15_min","Diss_20_min","Diss_25_min")])
TEST_JoBS <- cbind(ex_data_JoBS[ which(ex_data_JoBS$Group=='TEST'), ]
                   [c("Diss_15_min","Diss_20_min","Diss_25_min")])
equivalence_margin_EDNE_JoBS <- 297
test_EDNE_JoBS <- EDNE.EQ(X=REF_JoBS
                          , Y=TEST_JoBS
                          , eq_margin=equivalence_margin_EDNE_JoBS
                          , print.results = TRUE)


# Apart from simulation errors, a recalculation of the EDNE results 
# of some parts (normal distribution only) of the simulation study in 
# Hoffelder et al. (2015) can be done with the following code. Please note that
# the simulation takes approximately 7 minutes for 50.000 simulation 
# runs (number_of_simu_runs <- 50000). To shorten calculation time for 
# test users, number_of_simu_runs is set to 100 here and can/should be adapted.

library(MASS)
number_of_simu_runs <- 100
set.seed(2020)

mu1 <- c(41,76,97)
mu2 <- mu1 - c(10,10,10)
SIGMA_1 <- matrix(data = c(537.4 , 323.8 , 91.8 ,
                           323.8 , 207.5 , 61.7 ,
                           91.8  , 61.7  , 26.1) ,ncol = 3)
SIGMA_2 <- matrix(data = c(324.1 , 233.6 , 24.5 ,
                           233.6 , 263.5 , 61.4 ,
                           24.5  , 61.4  , 32.5) ,ncol = 3)
SIGMA   <- matrix(data = c(430.7 , 278.7 , 58.1 ,
                           278.7 , 235.5 , 61.6 ,
                           58.1  , 61.6  , 29.3) ,ncol = 3)


SIMULATION_SIZE_EDNE <- function(disttype , Hom , Var , mu_1 , mu_2 
                                  ,  n_per_group , n_simus ) {
  
  n_success_EDNE <- 0
  if 	( Hom == "Yes" ) 	{
    COVMAT_1 <- SIGMA
    COVMAT_2 <- SIGMA
  }
  else 	   {
    COVMAT_1 <- SIGMA_1
    COVMAT_2 <- SIGMA_2
  } 
  if 	( Var == "Low" ) 	{
    COVMAT_1 <- COVMAT_1 / 4
    COVMAT_2 <- COVMAT_2 / 4
  }
  
  d <- ncol(COVMAT_1)
  Mean_diff <- mu_1 - mu_2          # Difference of both exp. values
  dist_edne <- crossprod(Mean_diff) # true EDNE distance and equivalence margin
  
  if 	( n_per_group == 10 ) 	{
    cat("Expected value sample 1:",mu_1,"\n",
        "Expected value sample 2:",mu_2,"\n",
        "Covariance matrix sample 1:",COVMAT_1,"\n",
        "Covariance matrix sample 2:",COVMAT_2,"\n",
        "EM_EDNE:",dist_edne,"\n")
  }  
  
  for (i in 1:n_simus)  {
    if 	( disttype == "Normal" ) 	{
      REF <- mvrnorm(n = n_per_group, mu=mu_1, Sigma=COVMAT_1)
      TEST<- mvrnorm(n = n_per_group, mu=mu_2, Sigma=COVMAT_2)
    }
    n_success_EDNE  <- n_success_EDNE + EDNE.EQ.dissolution.profiles(X=REF
                                      , Y=TEST
                                      , print.results = FALSE)$testresult.num
  }
  empirical_succ_prob_EDNE  <- n_success_EDNE / n_simus  
  simuresults  <- data.frame(dist = disttype , Hom = Hom , Var = Var 
                            , dimension = d , em_edne = dist_edne   
                            , sample.size = n_per_group 
                            , empirical.size.edne = empirical_succ_prob_EDNE)
}

SIMULATION_LOOP_SAMPLE_SIZE <- function(disttype , Hom , Var , mu_1 , mu_2 
                                        ,  n_simus ) {
  
  run_10  <- SIMULATION_SIZE_EDNE(disttype = disttype , Hom = Hom , Var = Var 
                                  , mu_1 = mu_1 , mu_2 = mu_2 
                                  , n_per_group = 10 , n_simus = n_simus) 
  run_30  <- SIMULATION_SIZE_EDNE(disttype = disttype , Hom = Hom , Var = Var 
                                  , mu_1 = mu_1 , mu_2 = mu_2 
                                  , n_per_group = 30 , n_simus = n_simus) 
  run_50  <- SIMULATION_SIZE_EDNE(disttype = disttype , Hom = Hom , Var = Var 
                                  , mu_1 = mu_1 , mu_2 = mu_2 
                                  , n_per_group = 50 , n_simus = n_simus) 
  run_100 <- SIMULATION_SIZE_EDNE(disttype = disttype , Hom = Hom , Var = Var 
                                  , mu_1 = mu_1 , mu_2 = mu_2 
                                  , n_per_group = 100 , n_simus = n_simus) 
  RESULT_MATRIX <- rbind(run_10 , run_30 , run_50 , run_100)
  RESULT_MATRIX                               
}

simu_1 <- SIMULATION_LOOP_SAMPLE_SIZE(disttype = "Normal", Hom = "Yes" 
                                      , Var = "High" 
                                      , mu_1 = mu1  
                                      , mu_2 = mu2 
                                      , n_simus = number_of_simu_runs)
simu_2 <- SIMULATION_LOOP_SAMPLE_SIZE(disttype = "Normal", Hom = "Yes" 
                                      , Var = "Low" 
                                      , mu_1 = mu1  
                                      , mu_2 = mu2 
                                      , n_simus = number_of_simu_runs)
simu_3 <- SIMULATION_LOOP_SAMPLE_SIZE(disttype = "Normal", Hom = "No" 
                                      , Var = "High" 
                                      , mu_1 = mu1  
                                      , mu_2 = mu2 
                                      , n_simus = number_of_simu_runs)
simu_4 <- SIMULATION_LOOP_SAMPLE_SIZE(disttype = "Normal", Hom = "No" 
                                      , Var = "Low" 
                                      , mu_1 = mu1  
                                      , mu_2 = mu2 
                                      , n_simus = number_of_simu_runs)

FINAL_RESULT <- rbind(simu_1 , simu_2 , simu_3 , simu_4)

cat("****  Simu results n_simu_runs: ",number_of_simu_runs,"   **** \n")
FINAL_RESULT

}

