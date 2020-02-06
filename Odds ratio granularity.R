# Libraries
library("RMySQL")
library(tidyverse)
library(networkD3)
library(igraph)
library(plotly)

# Functions
dbDisconnectAll <- function(){
  ile <- length(dbListConnections(MySQL())  )
  lapply( dbListConnections(MySQL()), function(x) dbDisconnect(x) )
  cat(sprintf("%s connection(s) closed.\n", ile))
}

# General SQL query
query <- function(q)
{
  # Remote DB with password
  con <- dbConnect(MySQL(), 
                   user  = "guest",
                   password    = "guest",
                   dbname="meadows",
                   port = 3306,
                   host   = "sxouse.ddns.net")
  rs1 = dbSendQuery(con, q)
  return(as_tibble(fetch(rs1, n=-1)))
  dbDisconnectAll()
}

# Load the database
GetTheData <-  function()
{
  # GET DATA FROM DB
  # Remote DB with password
  con <- dbConnect(MySQL(), 
                   user  = "guest",
                   password    = "guest",
                   dbname="meadows",
                   port = 3306,
                   host   = "sxouse.ddns.net")
  
  
  q <- sprintf('select assembly_id, assembly_name, quadrat_count, community, quadrat_id, quadrat_size, visit_date, records_id, species.species_id, 
    species.species_name from assemblies
      join quadrats on quadrats.assembly_id = assemblies_id
      join visit_dates on quadrats.vd_id = visit_dates.vds_id
      join records on records.quadrat_id = quadrats_id
      join species on species.species_id = records.species_id
      # Two assemblies have 0 quadrat count; exclude A.capillaris_stolonifera;
      # exclude some odd assemblies with no assigned community
    where quadrat_count = 5 and species.species_id != 4 and community is not null
    and quadrat_size = "2x2";') 
  # NOTE: this extract includes "MG5", i.e. some MG5 communities where 
  # the team have not decided
  # on a sub-group.
  
  rs1 = dbSendQuery(con, q)
  return(as_tibble(fetch(rs1, n=-1)))
  dbDisconnectAll()
}

OddsRatio <-  function(d, A, B) # quadrat/assembly data, species_name, species_name
{
  d <- d %>% rename("id" = 1)
  A <- A[1]
  B <- B[1]
  t <- d %>% filter(species_name %in% c(A, B))
  As <- t %>% filter(species_name == A)
  Bs <- t %>% filter(species_name == B)
  j1 <- full_join(As, Bs, by = "id")
  # Get all the assemblies/quadrats, including ones with neither A nor B
  q <- d %>% distinct(id) 
  j2 <- (left_join(q, j1, by = "id") 
         # NOTE: column length q >= column length j1
         %>% mutate(AandB = !is.na(species_name.x) & !is.na(species_name.y))
         %>% mutate(AnotB = !is.na(species_name.x) & is.na(species_name.y))
         %>% mutate(BnotA = !is.na(species_name.y) & is.na(species_name.x))
         %>% mutate(neither = is.na(species_name.x) & is.na(species_name.y)))
  s <- colSums(j2[,4:7])
  # Standard Error https://en.wikipedia.org/wiki/Odds_ratio
  se <- sqrt((1/s[1] + (1/s[2]) + (1/s[3]) + (1/s[4]))) # std error of the log(o_r)
  o_r <- log((s[1]*s[4])/(s[2]*s[3])) # Return the LOG OR
  ci_low <- o_r - 1.96*se
  ci_high <- o_r + 1.96*se
  # Marginal sums of the joint probability distribution
  marg_x <- (s[1] + s[2])/(s[1] + s[2] + s[3] + s[4])
  marg_y <- (s[1] + s[3])/(s[1] + s[2] + s[3] + s[4])
  retval <- c(o_r, ci_low, ci_high, marg_x, marg_y)
  names(retval) <- c("odds_ratio", "ci_low", "ci_high", "marg_x", "marg_y")
  return(retval)
}

QuadratORGivenAssemblyOR <- function(r, marg_x, marg_y, sim_length = 10000, pool_size = 5)
  # r: assembly log odds ratio
  # marg_x: assembly marginal probability for species X
  # marg_y: assembly marginal probability for species Y
  # sim_length: number of simulated trials
  # pool_size: quadrat count per assembly
{
  # https://en.wikipedia.org/wiki/Odds_ratio
  # The odds ratio is a function of the cell probabilities, and conversely, 
  # the cell probabilities can be recovered given knowledge of the odds 
  # ratio and the marginal probabilities ...
  # Don't really need this as we could just retain the colSums from the
  # Odds ratio function and get the joint probability table directly
  
  r <- exp(r)
  s <- sqrt((1 + (marg_x + marg_y)*(r - 1))^2 + 4*r*(1 - r)*marg_x*marg_y)
  # Calculate the joint probability (contingency) table
  p11 <- ((marg_x + marg_y)*(r-1) - s + 1)/(2*(r-1))
  p10 <- marg_x - p11
  p01 <- marg_y - p11
  p00 <- 1 - marg_y - p10
  # Check we have actually recoverd
  # the assembly odds ratio, so the joint probability table is OK
  # r_check <- (p00*p11)/(p01*p10)
  # cat("check: ", r, r_check, "\n", sep = ", ")
  
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
  csums <- colSums(quads[5:8])
  quad_or <- log((csums[1]*csums[4])/(csums[2]*csums[3]))
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
  return (list(quad_or, ass_or))
} # end function


######
the_data <- GetTheData()
# Restrict analysis to species with more than 20 hits
sp_count <- the_data %>% group_by(species_id) %>% summarise(n = n())
sp_list <- sp_count %>% filter(n > 20)
# the_data with reduced species count
the_data <- the_data %>% filter(species_id %in% sp_list$species_id)
### the_data <- the_data %>% sample_n(1000)
rm(sp_count, sp_list)

# Select followed by distinct gives us just the occurrence of each species in each assembly,
# without multiple records of its occurrence.
# This stage was missing in the previous attempt at assembly_level odds ratio.
assembly_data <- the_data %>% select(assembly_id, species_name) %>% distinct()

quadrat_data <- the_data %>% select(quadrat_id, species_name)
# In order to use the odds ratio function, we need the species pairs
# It will be better to get the species pairs from the quadrat data: assemblies
# may have shared species that don't share a quadrat. All species sharing a 
# quadrat imply a shared assembly.
edges <- (full_join(quadrat_data, quadrat_data, by = "quadrat_id")
          %>% rename(from = species_name.x)
          %>% rename(to = species_name.y)
          %>%  select(-quadrat_id)
          %>% group_by(from, to) 
          %>% summarise(share_2x2 = n())
          %>% filter(from != to))
edges <- edges %>% filter(share_2x2 > 20)
nodes <- distinct(edges, from)
net <- graph_from_data_frame(d = edges, vertices = nodes, directed = F)
n2 <- simplify(net, edge.attr.comb = "first")
edges <- as_data_frame(n2, what="edges")
rm(nodes, net, n2)

# For each species pair, get the assembly odds ratio.
n <- length(row.names(edges))
assembly_pwor <- tibble(
  odds_ratio = 1:n,
  ci_low = 1:n,
  ci_high = 1:n,
  marg_x = 1:n,
  marg_y = 1:n)
for (i in seq_along(row.names(edges)))
{
  # assembly_pwor[i,1:5] <- AssemblyOddsRatio(edges$from[i], edges$to[i])
  assembly_pwor[i,1:5] <- OddsRatio(assembly_data, edges$from[i], edges$to[i])
} 

# Now take the same set of species pairs and make the quadrats odds ratios
quadrat_pwor <- tibble(
  odds_ratio = 1:n,
  ci_low = 1:n,
  ci_high = 1:n,
  marg_x = 1:n,
  marg_y = 1:n)
for (i in seq_along(row.names(edges)))
{
  quadrat_pwor[i,1:5] <- OddsRatio(quadrat_data, edges$from[i], edges$to[i])
}
quadrat_pwor <- bind_cols(edges, quadrat_pwor)
pwor <- bind_cols(quadrat_pwor, assembly_pwor)
colnames(pwor) <- c("from", "to", "share_2x2", "quadrat_or", "quadrat_ci_low", "quadrat_ci_high", 
                    "qmx", "qmy", "assembly_or", "assembly_ci_low", "assembly_ci_high", "amx", "amy")
rm(assembly_pwor, quadrat_pwor)
pwor <- pwor %>% filter(!is.infinite(assembly_or)) # %>% filter(assembly_or > 0.0)

# tibble to hold simulation results
expected_qor <- tibble(
  A = pwor$from,
  B = pwor$to,
  ci_low = pwor$quadrat_ci_low,
  ci_high = pwor$quadrat_ci_high,
  share_2x2 = pwor$share_2x2,
  obs = pwor$quadrat_or,
  expected_q = 0.0,
  expected_a = 0.0)
# Simulate quadrat OR expected if there were no local effect
for (i in seq_along(row.names(expected_qor)))
{
  simulated_or <- try(QuadratORGivenAssemblyOR(pwor$assembly_or[i], 
                      pwor$amx[i], pwor$amy[i], sim_length = 2000)) #2000))
  # Don't understand why as.numeric needed here.
  expected_qor[i,7] <- as.numeric(simulated_or[1])
  expected_qor[i,8] <- as.numeric(simulated_or[2])
}
expected_qor <- expected_qor %>% filter(!is.na(expected_q))

# Plot simulated vs observed qor
plt4 <- ggplot(expected_qor, aes(x=obs, y=expected_q)) +
  geom_abline(colour = "green")+
  geom_smooth(method = "auto", size = 0.5) +
  geom_point(aes(colour = share_2x2,text = paste(A, B, sep=",")), alpha = 0.5) +
  scale_colour_gradient(low = "sienna1", high = "black") +
  labs(x = "log(observed quadratOR)", y = "log(simulated quadratOR|assemblyOR)") +
  theme_grey() + coord_cartesian(xlim = c(-2, 3), ylim = c(-2,3))
plotly::ggplotly(plt4)

# write.csv(expected_qor, "expected QOR 10K iterations.csv", row.names = FALSE)
# 

# Check recovery of assembly odds ratio
# Plot observed aor vs aor recovered from simulation: should be straight line unit slope
recovered_aor <- bind_cols(pwor %>% select(assembly_or), expected_qor %>% select(expected_a))
# Don't understand why the as.numeric should be required here.
plt6 <- ggplot(recovered_aor, aes(x=assembly_or, y=as.numeric(expected_a)), colour = "blue") +
  geom_abline(colour = "blue") +
  geom_point() +
  labs(x = "log(observed assembly_or)", y = "log(recovered assembly_or)") +
  theme_grey() + coord_cartesian(xlim = c(-3, 3), ylim = c(-3,3))
plotly::ggplotly(plt6)

