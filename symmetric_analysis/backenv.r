require(ggplot2)

sc = 10000

b=0.1
s=200
be = rnorm(n=sc,b,sqrt(b/s))

rg= function(v) rgamma(sc,1/v,1/v)
#rg= function(v) rgamma(sc,(1+sqrt(4*v+1))/2/v+1,(1+sqrt(4*v+1))/2/v)
rexp = rexp(sc,log(2))
#rLN = function(v) rlnorm(sc,meanlog = -log(v+1)/2,sdlog =sqrt(log(v+1)))
rLN = function(v) rlnorm(sc,meanlog = 0,sdlog =sqrt(log(1/2 + 1/2 *sqrt(1 + 4 *v))))
rCN = function(v,a=0.1) {x=(1:60)/40;p=x[which.min(abs(c((exp(2*x)-exp(-2*x))/4/x-(exp(x)-exp(-x))^2/4/x^2)-v))];exp(runif(sc,-p,p))}
summary(rCN(1)); var(rCN(1));
summary(rLN(1)); var(rLN(1));
summary(rg(1))

lsd = function(be,b,r) {be=be[r>0.1]; r = r[r>0.1]; b/(be*r)-1}
logD = function(be,b,r) {log(b/(be*r))}

qplot(n, 
      data=rbind(
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
                 data.frame(m="LSD",r="Exponential",n=lsd(be,b,rexp),v="Var: 2.1"),
                 data.frame(m="LogD",r="Exponential",n=logD(be,b,rexp),v="Var: 2.1")
                 ),
      geom="density",color=m)+theme_classic()+facet_grid(v~r,scales="free")+ coord_cartesian(xlim=c(-4,4))+
  scale_color_brewer(name="",palette = "Dark2")+xlab("Penalty")+
  geom_vline(xintercept = 0,linetype=3)+theme(legend.position = c(0.87,0.2))

ggsave("backenv.pdf",width=10,height = 6)

