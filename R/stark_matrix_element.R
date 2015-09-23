#Calculates the matrix element for the Stark effect contribution to the energy of each state. Accepts an arbitrary n1, n2, l1, l2, j1, j2, mj1, mj2. Uses the adjusted Quantum Defect calculation.

stark_matrix_elem <- function(n1, n2, l1, l2, j1, j2, mj1, mj2){

  #Determines if the two mj terms are the same, if they are not then the matrix element is set to zero
  if(mj1 != mj2){
    StarkElem <- 0
  } else{
    #Determines if l2 = l1 +/-1. If not, the matrix element is set to zero.
    if((l2==(l1+1))|(l2==(l1-1))){

      #Calculates the spherical harmonic matrix element for mj1+1/2. If it is equal to zero, sets that term in the summation equal to zero.
      if(SphereMatElement(l1,l2,mj1+1/2)==0){

        SumPlus <- 0

      } else{

        SumPlus <- Clebsch_Gordan(l1,1/2,mj1+1/2,-1/2,j1,mj1)*Clebsch_Gordan(l2,1/2,mj1+1/2,-1/2,j2,mj1)*SphereMatElement(l1,l2,mj1+1/2)

      }
      #Calculates the spherical harmonic matrix element for mj1-1/2. If it is equal to zero, sets that term in the summation equal to zero.
      if(SphereMatElement(l1,l2,mj1-1/2)==0){

        SumMinus <- 0

      } else{

        SumMinus <- Clebsch_Gordan(l1,1/2,mj1-1/2,1/2,j1,mj1)*Clebsch_Gordan(l2,1/2,mj1-1/2,1/2,j2,mj1)*SphereMatElement(l1,l2,mj1-1/2)

      }

      #Calculates the stark matrix element.
      StarkElem <- RadialMatrixElement(1,n1,n2,l1,l2,j1,j2)*(SumPlus + SumMinus)

    } else{
      StarkElem <- 0
    }
  }

  StarkElem

}
