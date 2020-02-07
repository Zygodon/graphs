# Predict assembly odds ratio from quadrat odds ratio treating species hits combinations
# (X1Y1, X1Y0, X0Y1, X0Y0) as Bernouilli events and so calculating a 5-quadrat assembly
# joint probability matrix from the probability of at least one event on 5 trials: 
# 1 - (1-p)^5 giving a prediction of assembly level odds ratio based only on the quadrat 
# level associations.
# If quadrat level associations were unimportant we would expect 0 correlation.
# If assembly level correlations were unimportant we would expect perfect correlation.
# The plots suggest a strong contribution to assembly correlation from quadrat
# correlation.

# Libraries
library(tidyverse)
library(ggpmisc)

# Functions

p1in5 <- function(p, n=5)
{
  return(1 - (1-p)^n)
}

AssemblyORGivenQuadratOR <- function(jc)
# joint probability matrix                                     
{
  jp <- jc/sum(jc)
  a <- rep(0.0, 4)
  for (i in 1:4) a[i] <- p1in5(jp[i])
  lor <- log((a[1]*a[4])/(a[2]*a[3]))
  return(lor)
}

######

pwor <- read.csv("contingency.csv", header = TRUE, sep = ",")

# Make tibble for estimated assemblyOR from quadrat contingency
a <- pwor %>% select(A, B, share_2x2, aor, sfx)
a$aor_est <- 0.0
for (i in seq_along(row.names(a)))
{
  a$aor_est[i] <- AssemblyORGivenQuadratOR(c(pwor$jc1.x[i], pwor$jc2.x[i], pwor$jc3.x[i], pwor$jc4.x[i]))
} 
lin_model <- y ~ x # For stat_poly_eq

plt1 <- ggplot(a, aes(x=aor, y=aor_est), colour = "blue") +
  geom_abline(colour = "green") +
  geom_smooth(method = lm, size = 0.5) +
  stat_poly_eq(formula = lin_model,
               eq.with.lhs = "italic(hat(y))~`=`~",
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) +         
  geom_point(aes(colour = sfx, text = paste(A, B, sep=",")), alpha = 0.5) +
  scale_colour_manual(values = c("grey27", "sienna3")) +
  labs(x = "observed assembly log(odds ratio)", y = "predicted assembly log(odds ratio)") +
  theme_grey() + coord_cartesian(xlim = c(-2, 3), ylim = c(-2,3))
# plotly::ggplotly(plt1)
plot(plt1)


