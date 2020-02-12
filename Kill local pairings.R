# Kill local pairings
#
# Attempt to visualise the importance
# of specifically quadrat level association by reconstructing assembly
# level associations by simulation where the quadrat level association
# is destroyed. 
# In this case, if quadrat level were unimportant, the
# simuation would not affect the reconstructed assembly level odds ratios,
# so the graph of simulated assembly OR vs observed assembly OR would
# have the points on or near the unit slope.
# Conversely, if assembly level associations were unimportant (so, only
# quadrat associations contribute to the odds ratios), the points on the graph 
# would lie on a horizontal line crossing the origin. The result is near
# the second of these two scenarios, demonstrating the importance of 
# associations between plants in the same quadrat, not just in the 
# same assembly.

# The second plot checks that the simulation has in fact recovered the 
# observed assembly level odds ratios, while destroying any extra 
# association at the quadrat level.

# NOTE: Log odds ratio used throughout. i.e. for "odds ratio" everywhere, read
# "log(odds_ratio)".

# Libraries
library(tidyverse)
library(plotly)
library(ggpmisc)
library(DT)

# Functions

QuadratORGivenAssemblyOR <- function(r, marg_x, marg_y, sim_length = 10000, pool_size = 5)
  # r: assembly log odds ratio
  # marg_x: assembly marginal probability for species X
  # marg_y: assembly marginal probability for species Y
  # sim_length: number of simulated trials
  # pool_size: quadrat count per assembly
{

  r <- exp(r)
  s <- sqrt((1 + (marg_x + marg_y)*(r - 1))^2 + 4*r*(1 - r)*marg_x*marg_y)
  # Calculate the joint probability (contingency) table
  p11 <- ((marg_x + marg_y)*(r-1) - s + 1)/(2*(r-1))
  p10 <- marg_x - p11
  p01 <- marg_y - p11
  p00 <- 1 - marg_y - p10

  # Simulate the observed assembly data
  sim <-  tibble(trials = sample(c("X1Y1", "X0Y1", "X1Y0", "X0Y0"), sim_length, 
                                 prob = c(p11, p10, p01, p00), rep=T))

  # Recover simulated marginal totals
  simXY <- (sim %>% mutate(X = (trials == "X1Y1" | trials == "X1Y0"))
            %>% mutate(Y = (trials == "X1Y1" | trials == "X0Y1")))
  simXY$trial <- as.numeric(rownames(simXY))
  simXY <- simXY %>% select(trial, X, Y)
  
  # Make a tibble to put the quadrat_level simulation in
  quads <- tibble(q=rep(seq(1:pool_size), sim_length))
  quads$ass <- c( 1, 1 + seq(1:(pool_size*sim_length)) %/% pool_size)[1:(pool_size*sim_length)] # 1,1,1,1,1;2,2,2,2,2; ...
  # Set X and Y FALSE throughout
  quads$X <- F
  quads$Y <- F
  # THE KEY STEPS
  # Un-pool the samples, 5 to an "assembly", 
  # and assign the hits at random between the 5 "quadrats".
  # If there is a hit for X in "assembly" i, set one of the 5 "quadrats" X to TRUE at random
  # Similarly for Y. This should retain the assembly level association, but
  # completely destroy the local dependencies.
  for (i in 0:(sim_length-1))
  {
    k <- sample(1:pool_size,1) # A number in c(1,2,3,4,5)
    if(simXY$X[i+1]) quads$X[5*i+k] <- TRUE
    k <- sample(1:pool_size,1)
    if(simXY$Y[i+1]) quads$Y[5*i+k] <- TRUE
  }
  # Record the contingencies given no local effects 
  quads <- quads %>% mutate(X1Y1 = X & Y)
  quads <- quads %>% mutate(X1Y0 = X & !Y)
  quads <- quads %>% mutate(X0Y1 = !X & Y)
  quads <- quads %>% mutate(X0Y0 = !X & !Y)
  # Calculate the expected quadrat odds ratios given the observed assembly
  # odds ratios but with no local effects
  s <- colSums(quads[5:8])
  quad_or <- log((s[1]*s[4])/(s[2]*s[3]))
  se <- sqrt((1/s[1] + (1/s[2]) + (1/s[3]) + (1/s[4]))) # std error of the log(o_r)
  q_ci_low <- quad_or - 1.96*se
  q_ci_high <- quad_or + 1.96*se
  
  # Calculate the assembly OR: have we retained it?
  ass_x <-  quads %>% group_by(ass) %>% summarise(max(X))
  ass_y <-  quads %>% group_by(ass) %>% summarise(max(Y))
  ass <- left_join(ass_x, ass_y, by = "ass")
  ass <- ass %>% rename(X = "max(X)") %>% rename(Y = "max(Y)")
  ass <- (ass %>% mutate(X1Y1 = (X == 1)&(Y == 1))
              %>% mutate(X1Y0 = (X == 1)&(Y == 0))
              %>% mutate(X0Y1 = (X == 0)&(Y == 1))
              %>% mutate(X0Y0 = (X == 0)&(Y == 0))
              %>% select(X1Y1, X1Y0, X0Y1, X0Y0))
  csums <- colSums(ass)
  ass_or <- log((csums[1]*csums[4])/(csums[2]*csums[3]))
  return (list(quad_or, ass_or, q_ci_low, q_ci_high))
} # end function


######
# 

#pwor <- read.csv("contingency.csv", header = TRUE, sep = ",")

pwor$sim_q = 0.0
pwor$sim_a = 0.0
pwor$sim_q_low = 0.0
pwor$sim_q_high = 0.0

# Simulate OR expected if there were no local effect
for (i in seq_along(row.names(pwor)))
{
  s <- sum(pwor[i, 10:13])
  # assembly marginal probabilities
  amx <- (pwor$jc1.y[i] + pwor$jc3.y[i])/s
  amy <- (pwor$jc1.y[i] + pwor$jc2.y[i])/s
  simulated_or <- try(QuadratORGivenAssemblyOR(pwor$aor[i], 
                      amx, amy, sim_length = 2000)) #2000))
  # Don't understand why as.numeric needed here.
  pwor$sim_q[i] <- as.numeric(simulated_or[1]) # Destroyed quadrat ORs
  pwor$sim_a[i] <- as.numeric(simulated_or[2]) # Reconstructed assembly ORs
  pwor$sim_q_low[i] <- as.numeric(simulated_or[3]) #
  pwor$sim_q_high[i] <- as.numeric(simulated_or[4])  #
  cat("pair ", i, "\n")
}
pwor <- pwor %>% filter(!is.na(sim_q))
write.csv(pwor, "pwor.csv", row.names = FALSE)

# Plot simulated vs observed qor : expect near zero slope if
# local influence in strength of association was important.
lin_model <- y ~ x # For stat_poly_eq
plt4 <- ggplot(pwor, aes(x = qor, y = sim_q)) +
  geom_abline(colour = "green")+
  geom_smooth(method = lm, size = 0.5) +
  stat_poly_eq(formula = lin_model,
               eq.with.lhs = "italic(hat(y))~`=`~",
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) +         
  geom_point(aes(colour = sfx, text = paste(A, B, sep=",")), alpha = 0.5) +
  scale_colour_manual(values = c("grey27", "sienna3")) +
  labs(x = "observed quadrat odds ratio)", y = "reconstructed quadrat odds ratio)") +
  ggtitle("Reconstruct quadrat OR from jumbled assembly OR") +
  theme_grey() + coord_cartesian(xlim = c(-2, 3), ylim = c(-2,3))
plot(plt4)

# Plot observed quadrat OR ~ assembly OR, superimpose reconstructed
# quadrat OR
plt5 <- ggplot(pwor, aes(x=aor, y=qor)) +
  geom_smooth(method = lm) +
  geom_errorbar(aes(ymin = ci_low.x, ymax = ci_high.x), 
                size = 0.1, width = 0.1, colour = "green", alpha = 0.5) +
  geom_errorbar(aes(ymin = sim_q_low, ymax = sim_q_high), 
                size = 0.1, width = 0.1, colour = "black", alpha = 0.5) +
  geom_point(aes(x = pwor$aor, y = pwor$sim_q,
             text = paste(pwor$A, pwor$B, sep=",")), shape = 1, 
             colour = "grey20", alpha = 0.5) +
  geom_point(aes(colour = share_2x2, text = paste(A, B, sep=",")), alpha = 0.5) +
  scale_colour_gradient(low = "sienna1", high = "black") +
  labs(x = "log(Odds Ratio), assemblies", y = "log(Odds Ratio), quadrats") +
  theme_grey() + coord_cartesian(xlim = c(min(pwor$aor), max(pwor$qor)), 
                                 ylim = c(min(pwor$aor), max(pwor$qor)))
plotly::ggplotly(plt5) # 

friends <- datatable(pwor %>% filter((sim_q_high < ci_low.x) & sfx == "yes")
            %>% select(A, B, qor, share_2x2), 
            caption = 'Plants found sharing a quadrat',
            colnames = c('log(odds ratio)' = 4, 'shared quadrats' = 5))
friends <- friends %>%  formatSignif('log(odds ratio)',2)

enemies <- datatable(pwor %>% filter((sim_q_low > ci_high.x) & (sfx == "yes"))
            %>% select(A, B, qor, share_2x2),
            caption = 'Plants that avoid sharing a quadrat',
            colnames = c('log(odds ratio)' = 4, 'shared quadrats' = 5))
enemies <- enemies %>% formatSignif('log(odds ratio)',2)


