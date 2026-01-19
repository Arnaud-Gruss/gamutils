#' @title Relative importance of predictors in the model
#'
#' @description This function evaluates the relative importance of predictors in a model. It implements the approach
#' developed by Thuiller et al. (2012), which involves comparing the predictions of a fitted model
#' (the "full model") with the predictions of models after randomisation of a given predictor, i.e. after random
#' permutation of the values of a given predictor within the dataset used for the fitted model ("random models").
#' One minus the Pearson’s correlation coefficient between the predictions of the full model and the predictions
#' of a random model gives an indication of the relative importance of a given predictor in explaining the
#' response variable of the full model. Note that the Pearson’s correlation coefficient can be negative.
#' Thuiller et al. (2012) consider that these cases represent an even greater importance of the permutated predictor
#' in explaining the response variable than with a correlation of 0.
#'
#' @param Model model of interest.
#' @param Data data frame used to fit the model of interest.
#' @param Nrep the number of times the randomisation process is repeated (default: 10).
#' @param r_seed the seed value (default: 1234).
#'
#' @return relImportance a vector providing mean relative importance, the lower bound of relative importance
#' (considering a 95% confidence interval), median relative importance, and the upper bound of relative importance
#' (considenting a 95% confidence interval).
#'
#' @references
#' Thuiller W, Lafourcade B, Araujo M. 2012. The presentation manual for BIOMOD. Grenoble, France:
#' Laboratoire d’Écologie Alpine, Université Joseph Fourier.
#'
#' @export
relative_importance <- function( Model, Data, Nrep = 10, r_seed = 1234 ) {

  ######## Get the class and family of the model of interest
  model_class <- class( Model )[1]
  model_family <- Model$family$family

  ######## Get the response variable of the model and a list of the model predictors
  Vars <- dimnames( attributes( terms( formula( Model ), data = Data ) )$factors )[[1]]
  pTemp <- as.vector( unlist( strsplit( Vars, '[(,)]' ) ) )
  PTempF <- str_trim( unique( pTemp ), side = 'both' )
  response <- pTemp[1]
  Predictors <- names( Data )[!is.na( match( names( Data ),PTempF[-1] ) )]

  ######## Use a randomisation process to assess the relative importance of the different variables in the model
  relImportance <- matrix( nrow = length( Predictors ), ncol = 4 )
  mPerfOutTemp <- numeric( length = Nrep )
  rownames( mPerfOut ) <- Predictors
  colnames( mPerfOut ) <- c( 'mean', 'lower2.5', 'median', 'upper97.5' )
  PredRef <- try( predict( Model, newdata = Data, type = 'response'), silent = FALSE )
  set.seed( r_seed )

  for( i in 1 : length( Predictors ) ) {
    tempData <- Data[is.na( match( names( Data ), Predictors[i] ) )]

    for( j in 1 : Nrep ) {

      sIndex <- sample( 1 : nrow( tempData ), nrow( tempData ), replace = T )
      newData <- cbind( tempData, Data[sIndex,Predictors[i]] )
      names( newData )[ncol( newData )] <- Predictors[i]
      tempPred <- try( predict( Model, newdata = newData, type = 'response'), silent = FALSE )
      mPerfOutTemp[j] <- 1 - cor( PredRef, tempPred )

    }

    relImportance[i,] <- c( mean( mPerfOutTemp ), quantile( mPerfOutTemp, c( 0.025, 0.5, 0.975 ) ) )

  }

  return( relImportance )

}
