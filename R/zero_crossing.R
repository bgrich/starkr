#' Zero Crossing.
#'
#' \code{zero_cross} calculates the zero crossings for a pair of vectors x
#' and y.
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

#' Two-body Resonance Calculator.
#'
#' \code{resonance} calculates the energy crossover of a pair of atoms in state
#' A and a pair state of atoms in state B and state C
#'
#' This function calculates the energy \eqn{(2*E_A - (E_B + E_C))}. The function
#' requires a data frame with all of the states and their field values such that
#' Frame$Field gives the electric field values and Frame$Ecm gives the energy in
#' inverse centimeters (cm^-1). The second argument, AFrame, is a separate data
#' frame that already has the frame for state A selected. The other two inputs
#' should be strings that correspond to a the state column in Frame.
#'
#' @param Frame A data frame. A data frame containing all of the in the Stark
#'   map. Needs columns Field, Ecm, and state.
#' @param AFrame. A data frame. Needs columns Field and Ecm.
#' @param StateB. A character vector identifying state B.
#' @param StateC. A character vector identifying state C.
resonance <- function(Frame, AFrame, StateB, StateC){

  zero_cross(AFrame$Field, (2*AFrame$Ecm - ((Frame%>%filter(state %in% StateB)%>%select(Field, Ecm))$Ecm + (Frame%>%filter(state%in%StateC)%>%select(Field,Ecm))$Ecm)))

}
