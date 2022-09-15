setRefClass("linreg", fields=list())
  
  
linreg<-function(formula,data){
  
    # Finding x and y
    x<-model.matrix(formula, data=data)
    y<-as.matrix(data[all.vars(formula)[1]]) #as.vector does not work
    
    # Finding regressions coefficients
    regressions_coefficients <- solve( t(x)%*%x) %*% t(x)%*%y
    
    # Finding the fitted values
    fitted_values <- x %*% regressions_coefficients
    
    # Finding the residuals
    the_residuals <- y - fitted_values
    
    # Finding the degrees of freedom
    n <- nrow(x)
    p <- ncol(x)
    the_degrees_of_freedom <- n-p
    
    # Finding the residual variance
    the_residual_variance <- t(the_residuals) %*% the_residuals / the_degrees_of_freedom
    
    # Finding the variance of the regression coefficients
    the_variance_of_the_regression_coefficients <- drop(the_residual_variance) * solve((t(x)%*%x))
    
    # Finding the t values for each coefficient
    t_value <- regressions_coefficients / sqrt(diag(the_variance_of_the_regression_coefficients))
    p_value <- pt(regressions_coefficients,the_degrees_of_freedom)
    
    return(t_value)
  
  } 

