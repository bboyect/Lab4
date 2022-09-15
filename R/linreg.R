setRefClass("linreg", fields=list())

linreg<-function(formula,data){
  x<-model.matrix(formula, data=data)
  y<-as.matrix(data[all.vars(formula)[1]]) #as.vector does not work
  print(typeof(y))
  
  regressions_coefficients <- solve( t(x)%%x) %% t(x)%*%y
  fitted_values <- x %*% regressions_coefficients
  the_residuals <- y - fitted_values
  
  
  return(the_residuals)
  
}