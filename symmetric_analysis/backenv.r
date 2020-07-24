require(ggplot2)

sc = 100000

b=0.1
s=200


###################### Product ################

be = rnorm(n=sc,b,sqrt(b/s))

lsd = function(be,b,r) {be=be; r = r; b/(be*r)-1}
logD = function(be,b,r) {be=be; r =r; log(b/(be*r))}

#rg= function(v) rgamma(sc,(1+sqrt(4*v+1))/2/v+1,(1+sqrt(4*v+1))/2/v)
#rLN = function(v) rlnorm(sc,meanlog = -log(v+1)/2,sdlog =sqrt(log(v+1)))
rg= function(v) rgamma(sc,1/v,1/v)
rexp = function() stats::rexp(sc,log(2))
rLN = function(v) rlnorm(sc,meanlog = 0,sdlog =sqrt(log(1/2 + 1/2 *sqrt(1 + 4 *v))))
rCN = function(v,a=0.1) {x=(1:60)/40;p=x[which.min(abs(c((exp(2*x)-exp(-2*x))/4/x-(exp(x)-exp(-x))^2/4/x^2)-v))];exp(runif(sc,-p,p))}
summary(rCN(1/2)); var(rCN(1/2));
summary(rLN(1/2)); var(rLN(1/2));
summary(rg(1/2)); var(rg(1/2))

d=rbind(
  #data.frame(m="LSD",r="Strict",n=lsd(be,b,1),v=0),
  #data.frame(m="LogD",r="Strict",n=logD(be,b,1),v=0),
  data.frame(m="LSD",r="Gamma",n=lsd(be,b,rg(1/log(2)^2)),v="Var: 2.1"),
  data.frame(m="LogD",r="Gamma",n=logD(be,b,rg(1/log(2)^2)),v="Var: 2.1"),
  data.frame(m="LSD",r="Gamma",n=lsd(be,b,rg(1/3)),v="Var: 0.33"),
  data.frame(m="LogD",r="Gamma",n=logD(be,b,rg(1/3)),v="Var: 0.33"),
  data.frame(m="LSD",r="Gamma",n=lsd(be,b,rg(1/1000)),v="Strict"),
  data.frame(m="LogD",r="Gamma",n=logD(be,b,rg(1/1000)),v="Strict"),
  data.frame(m="LSD",r="LogNormal",n=lsd(be,b,rLN(1/1000)),v="Strict"),
  data.frame(m="LogD",r="LogNormal",n=logD(be,b,rLN(1/1000)),v="Strict"),
  data.frame(m="LSD",r="LogNormal",n=lsd(be,b,rLN(1/log(2)^2)),v="Var: 2.1"),
  data.frame(m="LogD",r="LogNormal",n=logD(be,b,rLN(1/log(2)^2)),v="Var: 2.1"),                
  data.frame(m="LSD",r="LogNormal",n=lsd(be,b,rLN(1/3)),v="Var: 0.33"),
  data.frame(m="LogD",r="LogNormal",n=logD(be,b,rLN(1/3)),v="Var: 0.33"),
  data.frame(m="LSD",r="LogUnif",n=lsd(be,b,rCN(1/log(2)^2)),v="Var: 2.1"),
  data.frame(m="LogD",r="LogUnif",n=logD(be,b,rCN(1/log(2)^2)),v="Var: 2.1"),  
  data.frame(m="LSD",r="LogUnif",n=lsd(be,b,rCN(1/3)),v="Var: 0.33"),
  data.frame(m="LogD",r="LogUnif",n=logD(be,b,rCN(1/3)),v="Var: 0.33"),
  data.frame(m="LSD",r="LogUnif",n=lsd(be,b,rCN(1/1000)),v="Strict"),
  data.frame(m="LogD",r="LogUnif",n=logD(be,b,rCN(1/1000)),v="Strict") , 
  data.frame(m="LSD",r="Exponential",n=lsd(be,b,rexp()),v="Var: 2.1"),
  data.frame(m="LogD",r="Exponential",n=logD(be,b,rexp()),v="Var: 2.1")
)
ggplot(aes(x=n), data=d)+
      geom_density(aes(color=m))+geom_histogram(aes(fill=m,y=..density..),position = "identity",binwidth = 0.05,alpha=0.4)+theme_classic()+facet_grid(v~r,scales="free")+ 
  scale_color_brewer(name="",palette = "Dark2")+scale_fill_brewer(name="",palette = "Dark2")+xlab("Penalty")+xlim(-10,10)+
  geom_vline(xintercept = 0,linetype=3)+theme(legend.position = c(0.87,0.2))+coord_cartesian(xlim=c(-4,4))
ggsave("product.pdf",width=10,height = 6)


ggplot(aes(x=n,fill=m), data=d)+
  geom_histogram(position = "identity",binwidth = 0.02,alpha=0.5)+theme_classic()+facet_grid(v~r,scales="free")+ 
  scale_fill_brewer(name="",palette = "Dark2")+xlab("Penalty")+xlim(c(-10,10))+
  geom_vline(xintercept = 0,linetype=3)+theme(legend.position = c(0.87,0.2))+coord_cartesian(xlim=c(-1,1))
ggsave("hist.pdf",width=10,height = 6)

######################### Compound


bg = function(v) {bt = rgamma(sc,1/v,1/v)*b; rnorm(n=5*sc,bt,sqrt(bt/s))}
gexp = function() {bt = rexp(sc,log(2))*b ; rnorm(n=5*sc,bt,sqrt(bt/s))}
bLN = function(v) {bt = rlnorm(sc,meanlog = 0,sdlog =sqrt(log(1/2 + 1/2 *sqrt(1 + 4 *v))))*b; rnorm(n=5*sc,bt,sqrt(bt/s))}
bCN = function(v,a=0.1) {
  x=(1:60)/40;
  p=x[which.min(abs(c((exp(2*x)-exp(-2*x))/4/x-(exp(x)-exp(-x))^2/4/x^2)-v))];
  bt = exp(runif(sc,-p,p)); 
  rnorm(n=5*sc,bt,sqrt(bt/s))}
summary(b/bCN(1/2)); var(b/bCN(1/2));
summary(bLN(1/2)); var(bLN(1/2));
summary(bg(1/2)); var(bg(1/2))


lsd = function(be,b,v) {b/be(v)-1}
logD = function(be,b,v) {log(b/be(v))}

