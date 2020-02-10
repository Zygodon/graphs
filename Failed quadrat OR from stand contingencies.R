# Attempt to predict expected quadrat OR with no local effects by 
# calcualating the cell probabilities in the conjoint quadrat 
# level probability table for each species pair from the conjoint assembly
# (stand) level probability tables.
# Does not give a meaningful result because of non-linearities in the process.


BernPr <- function(p, n)
{
  return(1-(1-p)^(1/n))
}

pwor <- read.csv("contingency.csv", header = TRUE, sep = ",")

pwor <- pwor %>% mutate(s = jc1.y + jc2.y + jc3.y + jc4.y)
pwor <- pwor %>% mutate(p1 = jc1.y/s) %>% mutate(b1 = BernPr(p1, 5))
pwor <- pwor %>% mutate(p2 = jc2.y/s) %>% mutate(b2 = BernPr(p2, 5))
pwor <- pwor %>% mutate(p3 = jc3.y/s) %>% mutate(b3 = BernPr(p3, 5))
pwor <- pwor %>% mutate(p4 = jc4.y/s) %>% mutate(b4 = BernPr(p4, 5))
pwor <- pwor %>% select(-s, -p1, -p2, -p3, -p4)
pwor <- pwor %>% mutate(predicted_or = (b1*b4)/(b2*b3))

# Plot simulated vs observed qor : expect near zero slope if
# local influence in strength of association was important.
lin_model <- y ~ x # For stat_poly_eq
plt7 <- ggplot(pwor, aes(x = qor, y = predicted_or)) +
  geom_abline(colour = "green")+
  geom_smooth(method = lm, size = 0.5) +
  stat_poly_eq(formula = lin_model,
               eq.with.lhs = "italic(hat(y))~`=`~",
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) +         
  geom_point(aes(colour = sfx, text = paste(A, B, sep=",")), alpha = 0.5) +
  scale_colour_manual(values = c("grey27", "sienna3")) +
  labs(x = "observed quadrat odds ratio)", y = "reconstructed quadrat odds ratio)") +
  ggtitle("Predict quadrat OR from assembly OR") +
  theme_grey() + coord_cartesian(xlim = c(-2, 3), ylim = c(-2,3))
plot(plt7)
