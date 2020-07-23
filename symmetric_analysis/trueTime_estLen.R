setwd("/Users/uym2/Exploration/LogDate/MBE_revision_2")

d = read.table("trueTime_estLen.txt",header=T)

require(ggplot2)

colors = c("LSD" = "#1B9E77", "LogD" = "#D95F02")

ggplot(d[d$estLen > 1e-3,]) + 
  geom_density(aes(x=log(mu*trueTime/estLen),color="LogD")) + 
  geom_density(aes(x=mu*trueTime/estLen-1,color="LSD")) + 
  geom_vline(xintercept = 0) + coord_cartesian(xlim=c(-3,10)) +
  scale_color_manual(values = colors) + 
  xlab("Penalty") + 
  facet_wrap(~clock,scale="free") + theme_classic() +
  theme(legend.position = "bottom",legend.title = element_blank())

ggsave("trueTime_estLen.pdf",width=10,height=6)
