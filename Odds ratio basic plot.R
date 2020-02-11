
# Libraries
library("RMySQL")
library(tidyverse)
library(networkD3)
library(igraph)
library(plotly)
library(ggpmisc)

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

JointContingency <- function(d, A, B) # quadrat/assembly data, species_name, species_name
{
  d <- d %>% rename("id" = 1)
  A <- A[1]
  B <- B[1]
  As <- d %>% filter(species_name == A)
  Bs <- d %>% filter(species_name == B)
  j1 <- full_join(As, Bs, by = "id")
  # Get all the assemblies/quadrats, including ones with neither A nor B
  q <- d %>% distinct(id) 
  j2 <- (left_join(q, j1, by = "id") 
         # NOTE: column length q >= column length j1
         %>% mutate(X1Y1 = !is.na(species_name.x) & !is.na(species_name.y))
         %>% mutate(X1Y0 = !is.na(species_name.x) & is.na(species_name.y))
         %>% mutate(X0Y1 = !is.na(species_name.y) & is.na(species_name.x))
         %>% mutate(X0Y0 = is.na(species_name.x) & is.na(species_name.y)))
  s <- colSums(j2[,4:7])
  return(s)
}

OddsRatio <-  function(s) # Joint contingency (colSums)
{
  # Standard Error https://en.wikipedia.org/wiki/Odds_ratio
  se <- sqrt((1/s[1] + (1/s[2]) + (1/s[3]) + (1/s[4]))) # std error of the log(o_r)
  o_r <- log((s[1]*s[4])/(s[2]*s[3])) # returning LOG OR
  ci_low <- o_r - 1.96*se
  ci_high <- o_r + 1.96*se
  retval <- c(o_r, ci_low, ci_high)
  names(retval) <- c("odds_ratio", "ci_low", "ci_high")
  return(retval)
}

p1in5 <- function(p, n=5)
{
  return(1 - (1-p)^n)
}

AssemblyORGivenQuadratOR <- function(pr)
# joint probability matrix                                     
{
  a <- rep(0.0, 4)
  for (i in 1:4) a[i] <- p1in5(pr[i])
  lor <- log((a[1]*a[4])/(a[2]*a[3]))
  return(lor)
}

SfChi <-  function(jc)
{
  x <- matrix(unlist(jc), ncol = 2, nrow = 2, byrow = T)
  ifelse(chisq.test(x)$p.value < 0.05, "yes", "no")
}

######

if (file.exists("contingency.csv"))
{
  pwor <- read.csv("contingency.csv")
} else {
  the_data <- GetTheData()
  # Restrict analysis to species with more than 20 hits
  sp_count <- the_data %>% group_by(species_id) %>% summarise(n = n())
  sp_list <- sp_count %>% filter(n > 20)
  # the_data with reduced species count
  the_data <- the_data %>% filter(species_id %in% sp_list$species_id)
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
  assembly_pwor <- tibble(
    A = edges$from,
    B = edges$to,
    jc1 = 0,
    jc2 = 0,
    jc3 = 0,
    jc4 = 0,
    aor = 0,
    ci_low = 0,
    ci_high = 0)
  
  for (i in seq_along(row.names(edges)))
  {
    s <- JointContingency(assembly_data, edges$from[i], edges$to[i])
    assembly_pwor[i, 3:6] <- s[1:4]
    assembly_pwor[i, 7:9] <- OddsRatio(s)
  } # Don't remove NAs at this stage
  
  
  # Now take the same set of species pairs and make the quadrats odds ratios
  quadrat_pwor <- tibble(
    A = edges$from,
    B = edges$to,
    jc1 = 0,
    jc2 = 0,
    jc3 = 0,
    jc4 = 0,
    qor = 0,
    ci_low = 0,
    ci_high = 0)
  
  for (i in seq_along(row.names(edges)))
  {
    s <- JointContingency(quadrat_data, edges$from[i], edges$to[i])
    quadrat_pwor[i, 3:6] <- s[1:4]
    quadrat_pwor[i, 7:9] <- OddsRatio(s)
  } # Don't remove NAs at this stage
  
  # quadrat_pwor <- bind_cols(edges, quadrat_pwor)
  pwor <- (left_join(quadrat_pwor, assembly_pwor, by = c("A", "B"))
          %>% filter(!is.na(qor))
          %>% filter(!is.infinite(qor))
          %>% filter(!is.na(aor))
          %>% filter(!is.infinite(aor)))
  rm(assembly_pwor, quadrat_pwor)
  
  # Add Chisquare test result (p < 0.05)
  sfx <- tibble( sfx = pwor$A) # Arbitrary character string
  for (i in seq_along(row.names(pwor)))
  {
     x <- matrix(unlist(pwor[i, 3:6]), ncol = 2, nrow = 2, byrow = T)
     sf <- SfChi(x) #sfx: significance Chisquare.test: p < 0.05
     sfx$sfx[i] <- sf
     cat(i, sf, "\n")
  }
  pwor$sfx <- sfx$sfx
  rm(sfx)
  
  # Join /share_2x2 for later use
  pwor <- left_join(pwor, edges, by = c("A" = "from", "B" = "to"))
  write.csv(pwor, "contingency.csv", row.names = FALSE)
}
lin_model <- y ~ x # For stat_poly_eq

# Plot quadrat OR vs assembly OR. Note slope < 1.0; assemblies OR
# tend to be higher than quadrat OR because of increased expectation that
# a pair of plants will be found together in the larger area.
plt2 <- ggplot(pwor, aes(x=aor, y=qor)) +
  geom_smooth(method = lm) + 
  stat_poly_eq(formula = lin_model,
               eq.with.lhs = "italic(hat(y))~`=`~",
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) +         
  geom_point(aes(colour = sfx, text = paste(A, B, sep=",")), alpha = 0.5) +
  scale_colour_manual(values = c("grey27", "sienna3")) +
  labs(x = "log(Odds Ratio), assemblies", y = "log(Odds Ratio), quadrats") +
  theme_grey() + coord_cartesian(xlim = c(min(pwor$aor), max(pwor$qor)), 
                                 ylim = c(min(pwor$aor), max(pwor$qor)))
# plotly::ggplotly(plt2) # plotly not available with stat_poly_eq
plot(plt2)

