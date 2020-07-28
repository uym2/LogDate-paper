# phyML estimated

d = read.table("trueTime_estLen.txt",header=T)

require(ggplot2)

colors = c("LSD" = "#1B9E77", "LogD" = "#D95F02")

ggplot(d[d$estLen > 1e-3,]) + 
  geom_density(aes(x=mu*trueTime/estLen-1,fill="LSD",color="LSD"),alpha=0.4) + 
  geom_density(aes(x=log(mu*trueTime/estLen),fill="LogD",color="LogD"),alpha = 0.4) + 
  geom_vline(xintercept = 0,linetype=3) +
  scale_fill_manual(values = colors) + 
  scale_color_manual(values = colors) + 
  xlab("Penalty (without square)") + 
  xlim(-10,10)+coord_cartesian(xlim=c(-3.3,3.3))+
  facet_wrap(~clock,scale="free") + theme_classic() +
  theme(legend.position = "None")

ggsave("trueTime_estLen.pdf",width=4,height=4)

###################################################

# Normal model of the estimated branches
sc = 500000
b=0.1
s=500

bg = function(v) { ee( rgamma(sc,1/v,1/v)*b)}
bexp = function() {ee (stats::rexp(sc,1)*b )}
bLN = function(v) {ee (rlnorm(sc,meanlog = 0,sdlog =sqrt(log(1/2 + 1/2 *sqrt(1 + 4 *v))))*b)}
bCN = function(v,a=0.1) {
  x=(1:60)/40;
  p=x[which.min(abs(c((exp(2*x)-exp(-2*x))/4/x-(exp(x)-exp(-x))^2/4/x^2)-v))];
  ee( exp(runif(sc,-p,p))*b)}

bgl2 = bg(1/6);   
bLNl2 = bLN(1/6); 
bCNl2 = bCN(1/6);
bex = bexp();

ee  = function(bt) vapply(rnorm(n=t,bt,sqrt(bt/s)),FUN=function(x) max(x,0.00001),1) 
ee  = function(bt) {a=rnorm(n=t,bt,sqrt(bt/s));a[a>0]}

bs=ee(1*b)

lsd = function(be,b) {b/be-1}
logD = function(be,b) {log(b/be)}

d1=rbind(
  data.frame(m="LSD",r="Gamma",n=lsd(bgl2,b),v="Var: 1/6",o=bgl2),
  data.frame(m="LogDate",r="Gamma",n=logD(bgl2,b),v="Var: 1/6",o=bgl2),
  data.frame(m="LSD",r="LogNormal",n=lsd(bLNl2,b),v="Var: 1/6",o=bLNl2),
  data.frame(m="LogDate",r="LogNormal",n=logD(bLNl2,b),v="Var: 1/6",o=bLNl2),
  data.frame(m="LSD",r="Exponential",n=lsd(bex,b),v="Var: 1",o=bex),
  data.frame(m="LogDate",r="Exponential",n=logD(bex,b),v="Var: 1",o=bex),
  data.frame(m="LSD",r="Strict",n=lsd(bs,b),v="Var: 0",o=bs),
  data.frame(m="LogDate",r="Strict",n=logD(bs,b),v="Var: 0",o=bs)
);

ggplot(aes(x=n), data=d1)+
  geom_density(aes(color=m))+geom_histogram(aes(fill=m,y=..density..),position = "identity",binwidth = 0.05,alpha=0.4)+
  theme_classic()+facet_wrap(~r,scales="free")+ 
  scale_color_brewer(name="",palette = "Dark2")+
  scale_fill_brewer(name="",palette = "Dark2")+xlab("Penalty (without square)")+xlim(-10,10)+
  geom_vline(xintercept = 0,linetype=3)+
  theme(legend.position = "bottom")+
  coord_cartesian(xlim=c(-3.3,3.3))
ggsave("compound_4panels.pdf",width=4,height = 4)

