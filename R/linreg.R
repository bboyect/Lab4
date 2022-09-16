# Defining the class LinReg  
LinReg <- setRefClass("LinReg", fields = list(t_values = "matrix", p_values = "matrix", regressions_coefficients = "matrix"),
                      methods = list(
                      show = function(){
                        names <- rownames(regressions_coefficients)
                        name <- colnames(regressions_coefficients)
                        output <- list(regressions_coefficients) 
                        print("hey eric")
                      }  
))


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
    t_values <- regressions_coefficients / sqrt(diag(the_variance_of_the_regression_coefficients))
    p_values <- 2*pt(t_values,the_degrees_of_freedom, lower.tail = FALSE)
    #linreg_object <- LinReg(t_values = t_values, p_values = p_values, regressions_coefficients = regressions_coefficients)
    
    return(p_values)
} 



