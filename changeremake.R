# remake slide from Chang 2018

library(ggplot2)

inc <- c(25,23,20,17,14,37,29,18,11,5,29,26,19,16,10)
wq  <- rep(1:5,3)
scn <- factor(rep(c('a. No Vaccination', ' b. Current Coverage', 'c. Equal Dose Distribution'), each = 5))
scn <- factor(scn, c('a. No Vaccination', ' b. Current Coverage', 'c. Equal Dose Distribution'))
d <- data.frame(inc,wq,scn)

ggplot(d, aes(y=inc, x=factor(wq))) + geom_bar(stat='identity') + theme_bw() + facet_wrap(~scn) +
  xlab('Wealth Quintile') + ylab("Proportion of Cases (%)") +
  theme(legend.position = "none",
        axis.title.y = element_text(size = 16),
        axis.text.y  = element_text(size = 12),
        axis.title.x = element_text(size = 16),
        axis.text.x  = element_text(size = 12))
