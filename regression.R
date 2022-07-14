cat("\014")      # clean console
library(dplyr)
library(readxl)
library(fixest)
table <- read_excel("regression_table.xlsx")
table <- rename(table, worker_ID = W_1)
table <- rename(table, firm_ID = W_2)
table <- rename(table, loq_wage_ID = W_3)


reg = feols(table$loq_wage_ID ~ 1 |worker_ID + firm_ID , data = table)

