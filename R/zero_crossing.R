#Determines the location of zero crossings for two arbitrary vectors. The vector y is the vector whose zero crossing are being measured. The vector x is an associated position vector for y.
zero_cross <- function(x, y){
  #Initiates variables to hold the index and x position for the zero crossing.
  ZeroPos <- numeric()
  ZeroX <- numeric()

  vector_sign <- sign(y)

  zeros <- which(vector_sign == 0)

  difference_vector <- diff(vector_sign)

  positive_cross <- which(difference_vector == 2)
  negative_cross <- which(difference_vector == -2)

  cross <- c(positive_cross, negative_cross)

  for(i in 1:length(cross)){
    if(abs(y[cross[i]]) > abs(y[cross[i]+1])){
      ZeroPos <- c(ZeroPos, cross[i]+1)
    } else {
      ZeroPos <- c(ZeroPos, cross[i])
    }
  }

  ZeroPos <- c(ZeroPos, zeros)
  ZeroX <- x[ZeroPos]

  #If there are no zero crossings a message is printed. Otherwise, a matrix is sent out with the index and x position of each zero crossing.
  if(length(ZeroPos) == 0){
    print("This vector contains no zero crossings.")
  } else {
    OutputMatrix <- cbind(ZeroPos, ZeroX)
    colnames(OutputMatrix) <- c("Index", "X.position")
    OutputMatrix
  }
}
