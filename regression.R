cat("\014")      # clean console
library(xlsx)
library(dplyr)
library(readxl)
library(fixest)
table <- read_excel("regression_table.xlsx")
table <- rename(table, worker_ID = W_1)
table <- rename(table, `worker productivity` = W_2)
table <- rename(table, firm_ID = W_3)
table <- rename(table, `firm productivity` = W_4)
table <- rename(table, stable = W_5)
table <- rename(table, loq_wage = W_6)

table$stable <- (table$stable == 1)

reg = feols(table$loq_wage ~ 1 | worker_ID + firm_ID , data = table)

fixeffs <- fixef(reg,sorted = TRUE)
table[,c("worker FE", "firm FE")] <- NA
for (n in 1:length(fixeffs[["firm_ID"]])) {
  table[["firm FE"]][table$firm_ID == n] <- fixeffs[["firm_ID"]][n]
}
for (n in 1:length(fixeffs[["worker_ID"]])) {
  table[["worker FE"]][table$worker_ID == n] <- fixeffs[["worker_ID"]][n]
}

write.xlsx(fixeffs[["firm_ID"]], file="firmFE.xlsx")

# headline
cor(table$`worker productivity`, table$`firm productivity`, method = "pearson", use = "complete.obs")
cor(table$`worker FE`, table$`firm FE`, method = "pearson", use = "complete.obs")
cor(table$`worker FE`, table$`worker productivity`, method = "pearson", use = "complete.obs")
cor(table$`firm FE`, table$`firm productivity`, method = "pearson", use = "complete.obs")
# stable matches
cor(table$`worker productivity`[table$stable], table$`firm productivity`[table$stable], method = "pearson", use = "complete.obs")
cor(table$`worker FE`[table$stable], table$`firm FE`[table$stable], method = "pearson", use = "complete.obs")
cor(table$`worker FE`[table$stable], table$`worker productivity`[table$stable], method = "pearson", use = "complete.obs")
cor(table$`firm FE`[table$stable], table$`firm productivity`[table$stable], method = "pearson", use = "complete.obs")
# unstable matches
cor(table$`worker productivity`[!table$stable], table$`firm productivity`[!table$stable], method = "pearson", use = "complete.obs")
cor(table$`worker FE`[!table$stable], table$`firm FE`[!table$stable], method = "pearson", use = "complete.obs")
cor(table$`worker FE`[!table$stable], table$`worker productivity`[!table$stable], method = "pearson", use = "complete.obs")
cor(table$`firm FE`[!table$stable], table$`firm productivity`[!table$stable], method = "pearson", use = "complete.obs")
