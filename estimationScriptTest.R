
# next function shifts the template by center and reevaluates at arbitrary points 
register<- function(x, Y, parAlphaOLD, parCenterOLD, template, DFBump ){
print("running script with  x, Y, parAlphaOLD, parCenterOLD, and template")
print("make sure these are correct!")

shiftTemplate <- function(x, center, template) {
   splint(template$x + center, template$y, x)
  
}
M<- ncol(Y) 
for (I in 1:5) {
  cat("iteration", I, fill=TRUE )
  parCenter<- parAlpha <- rep(NA, M)
  for (j in 1:M) {
     cat( j, " ")
    # estimate Center and scale with  fixed template   
    out <- nls(Y[, j] ~ 
                 alpha*shiftTemplate(x,  center, template ),
               start = (list(center = parCenterOLD[j], alpha=parAlphaOLD[j])
               )
    )
    parCenter[j] <- coefficients(out)[1]
    parAlpha[j] <- coefficients(out)[2]
  }
  ##### update parameters  
  RMSE1<- sqrt(mean( (parCenterOLD - parCenter)^2))
  parCenterOLD <- parCenter
  parAlphaOLD <- parAlpha
  ##### put data on common locations by centers and  also rescale  
  registerData <- matrix(NA, length(templateGrid), M)
  for (j in 1:M) {
    xCenter <- x - parCenter[j]
    # only consider shifted values that stay within template support.
    ind <- (templateGrid < max(xCenter)) & (templateGrid > min(xCenter))
    registerData[ind, j] <- splint(xCenter,
                                   Y[, j]/ parAlpha[j], templateGrid[ind])
  }
  #####  update bump estimate by mean of registered data  
  bumpRaw <- rowMeans(registerData, na.rm = TRUE)
  bumpNew<- sreg(templateGrid, bumpRaw, df=DFBump )$fitted.values
  bumpNew<- bumpNew/ max(bumpNew )
  
  RMSE2<- sqrt(mean( (bumpNew - template$y)^2))
  cat(I, "center ", RMSE1, "template ", RMSE2, fill=TRUE)
  
  # update the template
  template <- list(x = templateGrid,
                   y =  bumpNew)
}
print(" All done!")

return( 
  list(template=template, parCenter=parCenter, parAlpha= parAlpha )
  )
}
