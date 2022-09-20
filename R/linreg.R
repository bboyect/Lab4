#'  Linear Regression
#' 
#'  Run Linear Regression 
#' 
#' @param formula The formula
#' @param data The data, iris
#' @usage linreg(Petal.Length~Species, data = iris)
#' @return Linear regression for Iris
#' @export 
#' @importFrom ggplot2 ggplot aes geom_point geom_line labs ylab xlab theme element_text

# Defining the class LinReg  
LinReg <- setRefClass("LinReg", fields = list(t_values = "matrix", p_values = "matrix", regressions_coefficients = "matrix", the_residuals = "matrix", fitted_values = "matrix", the_residual_variance = "matrix", the_variance_of_the_regression_coefficients = "numeric", standard_error = "numeric",  sigma1 = "matrix" , the_degrees_of_freedom = "numeric"),
                      methods = list(
                      show = function(){
                        output <- drop(regressions_coefficients) 
                        print(output)
                      },
                      
                      medians_of_resdiuals = function(){
                        first_species_x <- median(fitted_values[1:50])
                        second_species_x <- median(fitted_values[51:100])
                        third_species_x <- median(fitted_values[101:150])
                        
                        first_species_y <- rep(median(the_residuals[1:50]), 50)
                        second_species_y <- rep(median(the_residuals[51:100]), 50)
                        third_species_y <- rep(median(the_residuals[101:150]), 50)

                        medians_data <- data.frame(x = c(first_species_x,second_species_x,third_species_x) , y = c(first_species_y, second_species_y, third_species_y))

                        return( medians_data)
                      },
                      
                      medians_of_standerdized_residuals = function(){
                        first_species_x <- median(fitted_values[1:50])
                        second_species_x <- median(fitted_values[51:100])
                        third_species_x <- median(fitted_values[101:150])
                        
                        sqrt_standardized_residuals <- sqrt(abs(the_residuals / sd(the_residuals)))
                        first_species_y <- rep(median(sqrt_standardized_residuals[1:50]), 50)
                        second_species_y <- rep(median(sqrt_standardized_residuals[51:100]), 50)
                        third_species_y <- rep(median(sqrt_standardized_residuals[101:150]), 50)
                        
                        medians_data <- data.frame(x = c(first_species_x,second_species_x,third_species_x) , y = c(first_species_y, second_species_y, third_species_y))
                        
                        return( medians_data)
                      },
                      
                      plot = function(){
                        median_values_x <- medians_of_resdiuals()[1]
                        median_values_y <- medians_of_resdiuals()[2]
                        median_values_y2 <- medians_of_standerdized_residuals()[2]
                        sqrt_standardized_residuals <- sqrt(abs(the_residuals / sd(the_residuals)))
                        test <- sqrt(abs(the_residual_variance))
                        testdata <- data.frame(fitted_values, the_residuals, test, median_values_x, median_values_y, median_values_y2, sqrt_standardized_residuals)
                        
                        plot1 <- ggplot() + 
                           geom_point(data=testdata, mapping =  aes(x=fitted_values, y= the_residuals)) +
                           geom_line(data=testdata, mapping =  aes(x=fitted_values, y= unlist(median_values_y)), color = "red")
                        
                        plot1 <- plot1 +  labs(title = "Residuals vs Fitted") +
                           theme(plot.title =  element_text(hjust = 0.5)) +  ylab("Residuals") +  xlab("Fitted values")
                        
                        plot2 <-  ggplot() + 
                           geom_point(data=testdata, mapping =  aes(x=fitted_values, y= unlist(sqrt_standardized_residuals))) +
                           geom_line(data=testdata, mapping =  aes(x=fitted_values, y= unlist(median_values_y2)), color = "red")
                        
                        plot2 <- plot2 +  labs(title = "Scale-Location") +
                           theme(plot.title =  element_text(hjust = 0.5)) +  ylab( expression(paste(sqrt("Standerdized residuals")))) +  xlab("Fitted values")
                        
                        plot_output <- list(plot1, plot2)
                        
                        return(plot_output)
                      },
                      
                      
                      resid = function(){
                        return(the_residuals)
                      },
                      
                      pred = function(){
                        return(fitted_values)
                      },
                      
                      coef = function(){
                        return(drop(regressions_coefficients))
                      },
                      
                      summary = function(){
                        
                        summary_output <- data.frame(regressions_coefficients, standard_error, t_values, p_values)
                        colnames(summary_output) <- c("regressions_coefficients", "standard_error", "t_values", "p_values")
                        print(summary_output)
                        cat("estimate of   sigma",  sigma1, "\n")
                        cat("Degree of freedom", the_degrees_of_freedom)
                      }
                      
                      
))


linreg<-function(formula,data){
  
    # Finding x and y
    x <<- model.matrix(formula, data=data)
    y <<- as.matrix(data[all.vars(formula)[1]]) #as.vector does not work
    
    # Finding regressions coefficients
    regressions_coefficients <<- solve( t(x)%*%x) %*% t(x)%*%y
    
    # Finding the fitted values
    fitted_values <<- x %*% regressions_coefficients
    
    # Finding the residuals
    the_residuals <<- y - fitted_values
    
    # Finding the degrees of freedom
    n <<- nrow(x)
    p <<- ncol(x)
    the_degrees_of_freedom <<- n-p
    
    # Finding the residual variance
    the_residual_variance <<- t(the_residuals) %*% the_residuals / the_degrees_of_freedom
    
    # Finding the variance of the regression coefficients
    the_variance_of_the_regression_coefficients <<- diag(drop(the_residual_variance) * solve((t(x)%*%x))) #use diag to eliminate unecessory covariances
    
    standard_error <<- sqrt(the_variance_of_the_regression_coefficients)
     sigma1 <<- sqrt(the_residual_variance)
    
    # Finding the t values for each coefficient
    t_values <<- regressions_coefficients / sqrt(the_variance_of_the_regression_coefficients)
    p_values <<- 2*pt(t_values,the_degrees_of_freedom, lower.tail = FALSE)
    linreg_object <<- LinReg(t_values = t_values, p_values = p_values, regressions_coefficients = regressions_coefficients,the_residuals = the_residuals, fitted_values = fitted_values, the_residual_variance = the_residual_variance, the_variance_of_the_regression_coefficients = the_variance_of_the_regression_coefficients , standard_error = standard_error , sigma1 =  sigma1 , the_degrees_of_freedom = the_degrees_of_freedom)
    
    return(linreg_object)
} 

mod_object <- linreg(Petal.Length~Species, data = iris)
mod_object$coef()
