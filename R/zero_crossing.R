#' Zero Crossing
#'
#' \code{zero_cross} calculates the zero crossings for a pair of vectors x
#' and y
#'
#' This function calculates the zero crossing for vectors x and y. The vector y
#' is the dependent variable and the vector x is the indepdenent variable.
#'
#' @param x A numeric. The independent variable. An associated position or time
#' vector to y.
#' @param y A numeric. The dependent variable. Typically a signal or function
#' of the independent variable x.
zero_cross <- function(x, y){
  # Initiates variables to hold the index and x position for the zero crossing.
  ZeroPos <- numeric()
  ZeroX <- numeric()

  # Determines the sign of the values in y
  vector_sign <- sign(y)

  # Determines any actual zeros
  zeros <- which(vector_sign == 0)

  # Calculates the difference in the sign
  difference_vector <- diff(vector_sign)

  # Determines the position of zero crossings that do not include zero itself
  positive_cross <- which(difference_vector == 2)
  negative_cross <- which(difference_vector == -2)

  cross <- c(positive_cross, negative_cross)

  # Looks at the zero crossings and determines which one is closest to zero
  for(i in 1:length(cross)){
    if(abs(y[cross[i]]) > abs(y[cross[i]+1])){
      ZeroPos <- c(ZeroPos, cross[i]+1)
    } else {
      ZeroPos <- c(ZeroPos, cross[i])
    }
  }

  ZeroPos <- c(ZeroPos, zeros)
  ZeroX <- x[ZeroPos]

  # If there are no zero crossings a message is printed. Otherwise, a matrix is
  # sent out with the index and x position of each zero crossing.
  if(length(ZeroPos) == 0){
    print("This vector contains no zero crossings.")
  } else {
    OutputMatrix <- cbind(ZeroPos, ZeroX)
    colnames(OutputMatrix) <- c("Index", "X.position")
    OutputMatrix
  }
}
