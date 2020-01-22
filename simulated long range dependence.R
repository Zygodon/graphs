quadrat_or_given_assembly <- function(r, marg_x, marg_y, sim_length = 10000, pool_size = 5)
#  r: assembly odds ratio
# marg_x: marginal probability for species X
# marg_y: marginal probability for species Y
# sim_length: number of simulated trials
# pool_size: quadrat count per assembly
{
   # https://en.wikipedia.org/wiki/Odds_ratio
  s <- sqrt((1 + (marg_x + marg_y)*(r - 1))^2 + 4*r*(1 - r)*marg_x*marg_y)
  p11 <- ((marg_x + marg_y)*(r-1) - s + 1)/(2*(r-1))
  p10 <- marg_x - p11
  p01 <- marg_y - p11
  p00 <- 1 - marg_y - p10
  # r_check <- (p00*p11)/(p01*p10)
  # sim_length <-  2000
  
  sim <-  tibble(trials = sample(c("X1Y1", "X0Y1", "X1Y0", "X0Y0"), sim_length, 
                                 prob = c(p11, p10, p01, p00), rep=T))
  # For the purpose of this function we don't need the assembly level odds ratio
  # joint_probs <- sim %>% group_by(trials) %>% count()
  # ass_or <-  (joint_probs[1,2]*joint_probs[4,2])/(joint_probs[2,2]*joint_probs[3,2])
  
  # Un-pool the samples, 5 to an "assembly", and assign the hits at random between the 5 "quadrats".
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
  # If there is a hit for X in "assembly" i, set one of the 5 "quadrats" X to TRUE at random
  # Similarly for Y. No short-range association here.
  for (i in 0:(sim_length-1))
  {
    k <- sample(1:pool_size,1)
    if(simXY$X[i+1]) quads$X[5*i+k] <- TRUE
    k <- sample(1:pool_size,1)
    if(simXY$Y[i+1]) quads$Y[5*i+k] <- TRUE
  }
  quads <- quads %>% mutate(X1Y1 = X & Y)
  quads <- quads %>% mutate(X1Y0 = X & !Y)
  quads <- quads %>% mutate(X0Y1 = !X & Y)
  quads <- quads %>% mutate(X0Y0 = !X & !Y)
  csums <- colSums(quads[5:8])
  quad_or <- (csums[1]*csums[4])/(csums[2]*csums[3])
  return (quad_or)
} # end function


