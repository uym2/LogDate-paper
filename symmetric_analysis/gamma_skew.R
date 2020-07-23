alpha = seq(3.5,10,0.1)

skew_lsd = 4*sqrt(alpha-2)/(alpha-3)
skew_logdate = trigamma(alpha)/(digamma(alpha)^(3/2))

mean_lsd = alpha/(alpha-1)-1
mean_logdate = digamma(alpha)-log(alpha)

d = rbind(data.frame(method="LSD",alpha=alpha,skewness=skew_lsd,mean=mean_lsd),
          data.frame(method="LogDate",alpha=alpha,skewness=skew_logdate,mean=mean_logdate))

require(ggplot2)

ggplot(d,aes(x=1/alpha,y=skewness,color=method)) + 
  xlab("gamma variance") + 
  geom_line() + theme_bw()
