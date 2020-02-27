
Rho <- function(sq, sa)
{
  b <- sq[1]/sum(sq)
  p <- sa[1]/sum(sa)
  pp <- 1 - ((1 - b)^5)
  return(pp/p)
}


tbl_rho <- tibble(A = pwor$A, B = pwor$B, rho = 0.0)
for (i in seq_along(row.names(pwor)))
{
  sq <- pwor[i, 3:6]
  sa <- pwor[i, 10:13]
  tbl_rho[i, 3] <- Rho(sq, sa)
}

tbl_rho <- mutate(tbl_rho, sc_rho = scale(rho))
a <- ggplot(tbl_rho, aes(sc_rho)) +
  geom_histogram(binwidth = 0.5)
plot(a)

friends <- tbl_rho %>%  filter(sc_rho > 1.95) %>% select(-sc_rho)
datatable(friends, caption = 'Plants that prefer sharing a quadrat') %>%  formatRound('rho',2)

enemies <- tbl_rho %>%  filter(sc_rho < -1.95) %>% select(-sc_rho)
datatable(enemies, caption = 'Plants that dislike sharing a quadrat') %>%  formatRound('rho',2)
