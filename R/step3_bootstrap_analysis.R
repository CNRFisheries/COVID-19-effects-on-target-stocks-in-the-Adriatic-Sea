# title: bootstrap analysis in support to paper "COVID-19-effects-on-target-stocks-in-the-Adriatic-Sea"
# author: Stefano Guicciardi
# date:  Last compiled on 23 of march 2022
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # set working directory to source file location

dati <- read.csv("../results/CMSY_estimate_combined.csv", header = TRUE, sep = ",")

dati$stock <- factor(dati$stock)

specie <- unique(dati$stock)

ls <- length(specie)

n <- 10000 # bootstrap number

mediaBBmsy <- rep(0, ls)

LCI_BB <- rep(0, ls)

UCI_BB <- rep(0, ls)

diffBBmsy <- rep(0, n)

b1 <- rep(0, n)

b2 <- rep(0, n)

set.seed(17)

for (i in 1:ls) {
  
  ds2019 <- dati %>% filter(stock == specie[i] & year == 2019)
  
  ds2020 <- dati %>% filter(stock == specie[i] & year == 2020)
  
  for (j in 1:n) {
    
    diffBBmsy[j] <- median(sample(ds2020$BBmsy, dim(ds2020)[1], replace = TRUE)) - median(sample(ds2019$BBmsy, dim(ds2019)[1], replace = TRUE))
    
  }
  
  mediaBBmsy[i] <- sprintf(mean(diffBBmsy), fmt = '%#.4f')
  
  LCI_BB[i] <- sprintf(quantile(diffBBmsy, 0.025), fmt = '%#.4f')
  
  UCI_BB[i] <- sprintf(quantile(diffBBmsy, 0.975), fmt = '%#.4f')
  
}

tab_BB <- cbind.data.frame(specie, mediaBBmsy, LCI_BB, UCI_BB) 

tab_BB

write.csv(tab_BB, "../results/bootstrap_results.csv", row.names = FALSE)
