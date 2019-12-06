

library(data.table)


n_mothers <- 1000000
males     <- rbinom(n_mothers,1,.5)

d <- data.table(mother = 1:n_mothers)

d[, birth1 := rbinom(.N,1,0.5)]

birthsremain <- sum(d$birth1)
b <- 2
while(birthsremain != 0){
  d[[paste0('birth',b)]] <- rbinom(n_mothers,1,0.5) 
  d[[paste0('birth',b)]][d[[paste0('birth',b-1)]]==1 | 
                         is.na(d[[paste0('birth',b-1)]])] = NA
  birthsremain <- sum(d[[paste0('birth',b)]],na.rm = TRUE)
  b <- b + 1
}

dd <- na.omit(melt(d, id.vars    = 'mother', 
                   value.name    = 'male', 
                   variable.name = 'birth'))
dd[, birth := as.numeric(gsub('birth','',birth))]

# sex split
mean(dd$male)

#avg number of kids per mother
# 1/lambda? exponential? = 2, lets see
nrow(dd)/n_mothers
