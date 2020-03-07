
# Make quadrat odds ratio and join wih aor
qor <- tibble(
  A = edges$from,
  B = edges$to,
  jc1 = 0,
  jc2 = 0,
  jc3 = 0,
  jc4 = 0,
  or = 0,
  ci_low = 0,
  ci_high = 0)

for (i in seq_along(row.names(edges)))
{
  s <- JointContingency(quadrat_data, edges$from[i], edges$to[i])
  qor[i, 3:6] <- s[1:4]
  qor[i, 7:9] <- OddsRatio(s)
}

pwor <- (left_join(aor, qor, by = c("A", "B"))
         %>% filter(!is.na(or.y))
         %>% filter(!is.infinite(or.y))
         %>% filter(p_val<0.05))

# Predict binomial p (probability of co-occurence in 5 quadrats) from bernouilli b
# (probability of co-occurence in 1 quadrat)
pred_p <- tibble(A = pwor$A,
             B = pwor$B,
             p = 0,
             b = 0,
             pp = 0)

# Joint contingency tables abcd; ja for stands, jq for quadrats
ja <- pwor %>% select(jc1.x:jc4.x)
jq <- pwor %>% select(jc1.y:jc4.y)

for (i in seq_along(row.names(pred_p)))
{
  p <- ja[i,1]/sum(ja[i,])
  b <- jq[i,1]/sum(jq[i,])
  pp <- 1 - (1 - b)^5
  
  pred_p$p[i] <- unlist(p)
  pred_p$b[i] <- unlist(b)
  pred_p$pp[i] <- unlist(pp)
}

rm(ja, jq)

plt1 <- ggplot(pred_p, aes(p, pp)) +
  geom_point() +
  geom_abline(colour = "red") +
  xlim(0, 1) +
  ylim(0, 1) +
  xlab("observed P(co-occurrence)") +
  ylab("P(co-occurence) estimated from single trial")
plot(plt1)
rm