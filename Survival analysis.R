### Survival analysis - glmnet and Kaplan-Meier
### glmnet.R
# Packages
library(Biobase)
library(GEOquery)
require(glmnet)
require(survival)

# Input file
data=read.table("Genes_expression.csv", row.names=1, header=TRUE, sep=";", stringsAsFactors=FALSE)
data$time <- as.numeric(data$time)
data$status <- as.numeric(data$status)

# cross-val
x = data.matrix(data[,1:(ncol(data)-3)])
y = Surv(data$time, data$status)

cvfit = cv.glmnet(x=x, y=y, family='cox', nfolds=nrow(data), grouped = T)
cvfit
plot(cvfit)

# Set nlambda based on output
fit <- glmnet(x=x, y=y, family='cox', nlambda=1, lambda = cvfit$lambda.1se)

coefs <- coef(fit, s = cvfit$lambda.1se)
df <- data.frame(name = coefs@Dimnames[[1]][coefs@i + 1], coefficient = coefs@x)


### Kaplan-Meier
# Packages
library(survival)
library(ranger)
library(dplyr)
library(ggfortify)
library(survminer)

# Input file
data <- read.table("Gene_expression.csv", row.names=1, header=TRUE, sep=";", stringsAsFactors=FALSE)

data$time <- as.numeric(data$time)
md<- median(data$score)
data <- mutate(data, status = ifelse((score <  md), 0, 1), status = as.numeric(status)) 

fit <- survfit(Surv(time, status) ~ Risk, data = data)

ggsurvplot(fit, xlab= "Time in days",
           pval = TRUE, conf.int = F,
           data = data, 
           title ="Gene_name" ,
           risk.table.col = "strata",
           linetype = "strata", 
           dpi=600,
           palette = c( "blue" , "green", "red"),
           font.family = "Arial",
           ggtheme = theme_classic2(base_size = 15), 
           font.tickslab = c(15),
           text = element_text(size = 15))
