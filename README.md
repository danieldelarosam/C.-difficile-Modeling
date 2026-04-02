# *Clostridioides difficile* Transmission Model

This repository contains the code supporting the analyses presented in the manuscript:

"Clostridioides difficile Colonization Amplification despite Limited In-Hospital Transmission: A Modeling Study."

---

## Repository Structure

- **1_Model_Optimization_and_MonteCarlo.R**  
  R script to run the baseline compartmental model simulating *C. difficile* transmission and colonization amplification in a healthcare setting. Includes optimization of the transmission rate, Monte Carlo simulations, and estimation of key epidemiological metrics.

- **2_Univariate_sensitivity_analysis.R**  
  R script to conduct a univariate sensitivity analysis by varying each model parameter individually, recalibrating the transmission scaling parameter (δ), and recording the resulting  intrinsic reproduction number (Ri) and the colonization amplification index (Ai).

- **3_PRCC_Sensitivity_Analysis.R**  
  R script to perform partial rank correlation coefficient (PRCC) analyses to identify key drivers of variation in the Ri and the Ai.

- **4_Sensitivity_Analysis_Pi.R**  
  R script to perform sensitivity analysis on community-onset case definitions (π values). This script evaluates the impact of different thresholds for defining healthcare-associated *C. difficile* infections on Ai and Ri.

- **5_Interventions_Analysis.R**  
  R script to evaluate the impact of different intervention strategies (screening, isolation, and prophylaxis) on *C. difficile* transmission. The script runs simulations under varying screening coverage (efe_2 parameter) and estimates reductions in both CDI incidence and Ai.

- **6_NextGen_Matrix_Ri_Calculation.m**  
  MATLAB script to compute the intrinsic reproduction number for *C. difficile* using the next-generation matrix approach. The script derives the next-generation matrix symbolically, substitutes model parameters, and simplifies the final expression.
