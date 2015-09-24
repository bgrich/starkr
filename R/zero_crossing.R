#Determines the location of zero crossings for two arbitrary vectors. The vector y is the vector whose zero crossing are being measured. The vector x is an associated position vector for y.

zero_cross1 <- function(x, y){
  #Initiates variables to hold the index and x position for the zero crossing.
  ZeroPos <- numeric()
  ZeroX <- numeric()

  for(i in 1:(length(y)-1)){
    #If the value of y[i] is zero then it is kept
    if(sign(y[i]) == 0){
      ZeroPos <- c(ZeroPos, i)
      ZeroX <- c(ZeroX, x[i])
      next
    }
    #If the sign changes from i to i+1, then the value is kept.
    if(sign(y[i])!=sign(y[i+1])){
      if(abs(y[i]) > abs(y[i+1])){
        ZeroPos <- c(ZeroPos, i+1)
        ZeroX <- c(ZeroX, x[i+1])
      } else {
        ZeroPos <- c(ZeroPos, i)
        ZeroX <- c(ZeroX, x[i])
      }
    }
  }
  #If there are no zero crossings a message is printed. Otherwise, a matrix is sent out with the index and x position of each zero crossing.
  if(length(ZeroPos) == 0){
    print("This vector contains no zero crossings.")
  } else {
    OutputMatrix <- cbind(ZeroPos, ZeroX)
    colnames(OutputMatrix) <- c("Index", "X.position")
    OutputMatrix
  }
}



zero_cross2 <- function(x, y){
  #Initiates variables to hold the index and x position for the zero crossing.
  ZeroPos <- numeric()
  ZeroX <- numeric()

  vector_sign <- sign(y)

  zero_pos <- which(vector_sign == 0)

  difference_vector <- diff(vector_sign)

  positive_cross <- which(difference_vector == 2)
  negative_cross <- which(difference_vector == -2)

  ZeroPos <- c(zero_pos, positive_cross, negative_cross)
  ZeroX <- c(x[zero_pos], x[positive_cross], x[negative_cross])

  #If there are no zero crossings a message is printed. Otherwise, a matrix is sent out with the index and x position of each zero crossing.
  if(length(ZeroPos) == 0){
    print("This vector contains no zero crossings.")
  } else {
    OutputMatrix <- cbind(ZeroPos, ZeroX)
    colnames(OutputMatrix) <- c("Index", "X.position")
    OutputMatrix
  }
}
