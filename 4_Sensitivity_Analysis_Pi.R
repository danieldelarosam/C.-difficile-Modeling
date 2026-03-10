## -------------------------------
## 0. Required Packages 
## -------------------------------

if (!require(pacman)) install.packages("pacman")
library(pacman)
p_load(deSolve, ggplot2, dplyr)

# IMPORTANT: Run Script 1 first.
# This script assumes the following objects were generated:
# - initial_conditions
# - times
# - parameters_list (with num_simulations iterations)
# - comp_model (ODE function)

## -------------------------------
## 1. Sensitivity Analysis for 'pi' Values
## -------------------------------

### 1.1 Define 'pi' Values for Sensitivity Analysis
pi_values <- c(0.003347184, 0.003924284, 0.004570637, 0.004986150, 0.005609418)

sensitivity_simulations <- data.frame()

### 1.2 Run Simulations for Each 'pi' Value
# For each value of 'pi', optimize delta and run simulations.

for (pi_value in pi_values) {
  
  parameters_list$pi <- rep(pi_value, num_simulations)
  
  for (i in 1:num_simulations) {
    current_parameters <- lapply(parameters_list, function(x) x[i])
    
    delta_optimization <- function(delta) {
      current_parameters['delta'] <- delta
      result <- ode(y = initial_conditions, times = times, func = comp_model, parms = current_parameters)
      final_incidence <- tail(result[,"dDiagnosed_CDI"], 1)
      abs(639 - final_incidence)
    }
    
    optimized_delta <- optimize(f = delta_optimization, interval = c(0, 1), tol = 1e-10)
    current_parameters['delta'] <- optimized_delta$minimum
    optimized_output <- ode(y = initial_conditions, times = times, func = comp_model, parms = current_parameters)
    
    final_step <- optimized_output[nrow(optimized_output), ]
    
    colonization_amplification <- round((as.numeric(current_parameters[["psi_3"]]) *
                                           sum(final_step[c("E", "C", "K_1", "K_2")])) /
                                          (as.numeric(current_parameters[["lambda"]]) * 
                                             as.numeric(current_parameters[["ene"]])), 1)
    
    temp_results <- data.frame(
      simulation = i,
      pi_value = pi_value,
      dAdmission = final_step["dAdmission"],
      dDiagnosed_CDI = final_step["dDiagnosed_CDI"],
      E_equi = final_step["E"],
      C_equi = final_step["C"],
      K1_equi = final_step["K_1"],
      K2_equi = final_step["K_2"],
      Colonization_amplification = colonization_amplification
    )
    
    param_data <- as.numeric(unlist(current_parameters))
    names(param_data) <- names(current_parameters)
    temp_results <- cbind(temp_results, t(param_data))
    
    sensitivity_simulations <- rbind(sensitivity_simulations, temp_results)
  }
}

## -------------------------------
## 2. Summary and Visualization of Amplification Index
## -------------------------------

### 2.1 Summary Statistics for Colonization Amplification Index
summary_variable <- function(data, column, group_var = "pi_value") {
  data %>%
    group_by(across(all_of(group_var))) %>%
    summarise(
      median = median(.data[[column]], na.rm = TRUE),
      Q1 = quantile(.data[[column]], 0.25, na.rm = TRUE),
      Q3 = quantile(.data[[column]], 0.75, na.rm = TRUE)
    )
}

summary_variable(sensitivity_simulations, "Colonization_amplification")


### 2.2 Boxplot of Colonization Amplification Index
plot_boxplot <- function(data, column, title, ylab, x_labels) {
  ggplot(data, aes(x = factor(pi_value), y = .data[[column]])) +
    geom_boxplot(fill = "gray", color = "black") +
    labs(title = title, x = "Time to define healthcare associated infection", y = ylab) +
    scale_x_discrete(labels = x_labels) +
    theme_minimal(base_size = 12)
}

plot_boxplot(sensitivity_simulations, "Colonization_amplification",
             "Amplification Ratio of Colonization", 
             "Index (Colonized Discharged / Admitted)",
             c("1 day", "2 days", "3 days", "4 days", "5 days"))

## -------------------------------
## 3. Estimation of Reproduction Number (R)
## -------------------------------

### 3.1 Reorganize Parameters for R Calculation
parameters_list_updated <- list(
  d = sensitivity_simulations$delta,
  x = sensitivity_simulations$equis,
  a = sensitivity_simulations$alfa,
  e = sensitivity_simulations$epsilon,
  v = sensitivity_simulations$ve,
  l = sensitivity_simulations$lambda,
  n = sensitivity_simulations$ene,
  pi = sensitivity_simulations$pi_value,
  m = sensitivity_simulations$eme,
  z = sensitivity_simulations$zeta,
  p_1 = sensitivity_simulations$psi_1,
  p_2 = sensitivity_simulations$psi_2,
  p_3 = sensitivity_simulations$psi_3,
  p_4 = sensitivity_simulations$psi_4,
  h_1 = sensitivity_simulations$ache_1,
  h_2 = sensitivity_simulations$ache_2,
  g_1 = sensitivity_simulations$gamma_1,
  g_2 = sensitivity_simulations$gamma_2,
  f_2 = sensitivity_simulations$efe_2,
  f_1 = sensitivity_simulations$efe_1,
  s_1 = sensitivity_simulations$sigma_1,
  s_2 = sensitivity_simulations$sigma_2,
  index = sensitivity_simulations$Colonization_amplification
)

### 3.2 Calculate R Using Next-Generation Matrix
R_expr <- expression((d*l*(a + p_1*z)*(p_3*p_4^2*x + p_4^2*v*x + e*f_1*p_3*v + e*h_1*p_3*v + 
                                         f_1*h_1*p_3*x + e*p_3*p_4*v + f_1*h_1*v*x + f_1*p_3*p_4*x + 
                                         h_1*p_3*p_4*x + f_1*p_4*v*x + h_1*p_4*v*x - e*p_4^2*v*x -
                                         e*f_1*g_1*p_3*v - e*f_1*p_4*v*x - e*h_1*p_4*v*x -
                                         e*f_1*h_1*s_1*v*x))/(p_2*p_3*(a + p_1)*(f_1 + p_4)*(h_1 + p_4)*(p_3 + v)))

calculate_R <- function(i) {
  params <- lapply(parameters_list_updated, `[[`, i)
  eval(R_expr, envir = as.list(params))
}

parameters_list_updated <- as.data.frame(parameters_list_updated)
parameters_list_updated$R <- sapply(seq_len(nrow(parameters_list_updated)), calculate_R)
sensitivity_simulations$R <- parameters_list_updated$R

## -------------------------------
## 4. Summary and Visualization of Reproduction Number (R)
## -------------------------------

### 4.1 Summary Statistics for Reproduction Number (R)
summary_variable(sensitivity_simulations, "R")

### 4.2 Boxplot of Reproduction Number (R)
plot_boxplot(sensitivity_simulations, "R",
             "Intrinsic Reproduction Number (R)", 
             "R",
             c("1 day", "2 days", "3 days", "4 days", "5 days")) +
  geom_hline(yintercept = 1, linetype = "dashed")
