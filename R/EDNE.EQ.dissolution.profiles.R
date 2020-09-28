EDNE.EQ.dissolution.profiles <-
function( X , Y , alpha=0.05 , print.results=TRUE){
  
  if (print.results) {
    cat("\n")
    cat("*******************************************************************************"  , "\n\n")
    cat("********************* \t The EDNE-test for equivalence \t***********************" , "\n")
    cat("***********           for comparing dissolution profiles           ************"  , "\n")
    cat("******* 10% acceptance criterion analogous to similarity factor f2 \t*******"     , "\n\n")
    cat("*******************************************************************************"  , "\n\n")
  }
  
  m_X         <- nrow(X)                                                # sample size first sample
  n_Y         <- nrow(Y)                                                # sample size second sample
  N           <- m_X + n_Y                                              # total sample size pooled sample
  d           <- ncol(X)                                                # number of variables / dimension
  Mean_X      <- colMeans(X)                                            # Mean of first sample
  Mean_Y      <- colMeans(Y)                                            # Mean of second sample 
  Mean_diff   <- Mean_X - Mean_Y                                        # Difference of both sample mean vectors
  
  S_X         <- var(X)                                                 # empirical covariance matrix of first sample
  S_Y         <- var(Y)                                                 # empirical covariance matrix of second sample
  
  D_ems       <- matrix(10,d,1)                                         # Vector for calculating equivalence margin according to
  eq_margin   <- t(D_ems) %*% D_ems                                     # similarity factor f2 (see Suarez et al. (2020))
  
  Q_EDNE      <- crossprod(Mean_diff)                                   # estimated Euclidean distance
  
  MATRIX_HELP <-  S_X/m_X + S_Y/n_Y
  asympt_var  <-  4*t(Mean_diff) %*% MATRIX_HELP %*% Mean_diff          # estimated asymptotic variance
  
  numerator   <- Q_EDNE -  eq_margin                                    # numerator of EDNE test statistic
  denominator <- sqrt(asympt_var)                                       # denominator of EDNE test statistic
  
  quantile		<- qnorm( p=alpha , mean = 0, sd = 1 )                    # quantile of standard normal distribution
  teststat		<- numerator / denominator 	                              # teststatistic
  
  if (print.results) {
    cat(" Summary statistics:"                                                           , "\n\n",    
        "Sample size of first sample X:"                   ,"\t" ,m_X                    ,   "\n",    
        "Sample size of second sample Y:"                  ,"\t" ,n_Y                    ,   "\n",
        "Sample size of pooled sample (X,Y):"              ,"\t" ,N                      ,   "\n",
        "Dimension (Number of variables):"                 ,"\t" ,d                      ,   "\n",
        "Mean of first sample:"                        ,"\t\t\t" ,Mean_X                 ,   "\n",
        "Mean of second sample:"                         ,"\t\t" ,Mean_Y                 ,   "\n",
        "Difference of both means:"                      ,"\t\t" ,Mean_diff              ,   "\n",
        "Empirical covariance matrix S_X:"                 ,"\t" ,round(S_X,digits=2)    ,   "\n",
        "Empirical covariance matrix S_Y:"                 ,"\t" ,round(S_Y,digits=2)    ,   "\n",
        "Estimated Euclidean distance:"                ,"\t\t"   ,Q_EDNE                 ,   "\n",      
        "Equivalence margin:"                          ,"\t\t\t" ,eq_margin              ,   "\n",      
        "asymptotic variance:"                         ,"\t\t\t" ,asympt_var             ,   "\n",   
        "Significance level:"                          ,"\t\t\t" ,alpha                  ,   "\n"
    )
  }
  
  if 	( teststat < quantile ) 	{
    erg_text	<-	'Equivalence comparison successful'	                  # Test result text 
    erg       <- 1                                                      # Test result  
  }
  else 	   {
    erg_text	<-	'Equivalence comparison not successful'		            # Test result text 
    erg       <- 0                                                      # Test result  
  }
  
  p_value			<- pnorm( q=teststat , mean = 0, sd = 1 )                 # p-value of EDNE-test for equivalence  
  
  if (print.results) {
    cat(" Teststatistic:"                                        ,"\t\t\t" ,teststat ,   "\n",    
        "Quantile of stand. normal dist.:"                       ,"\t"     ,quantile ,   "\n\n",
        "Decision in favor (1) or against (0) similarity of dissolution profiles: "               ,erg      ,   "\n\n",
        "Test result:"                                                               ,   "\n\n",
        "\t p-value of the EDNE-test for equivalence: p ="                 ,p_value  ,   "\n\n",
        "\t\t"                                                             ,erg_text ,   "\n\n",
        "******************************************************************************"  ,   "\n\n"
    )
  }
  result.summary <- data.frame(p.value = p_value , testresult.num = erg , testresult.text = erg_text )
  result.summary	          # Return of Test results  
}
