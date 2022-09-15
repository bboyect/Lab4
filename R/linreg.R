
  setRefClass("linreg", fields=list())
  
  linreg<-function(formula,data){
    x<-model.matrix(formula, data=data)
    y<-as.matrix(data[all.vars(formula)[1]]) #as.vector does not work
    
    regressions_coefficients <- solve( t(x)%*%x) %*% t(x)%*%y
    fitted_values <- x %*% regressions_coefficients
    the_residuals <- y - fitted_values
    n <- nrow(x)
    p <- ncol(x)
    the_degrees_of_freedom <- n-p
    the_residual_variance <- t(the_residuals) %*% the_residuals / the_degrees_of_freedom
    the_variance_of_the_regression_coefficients <- drop(the_residual_variance) * solve((t(x)%*%x))
    
    return(the_variance_of_the_regression_coefficients)
  
  }
  