require(ggplot2)
require(reshape2)
require(scales)


sc = 500000

b=0.1
s=500


###################### Product ################

be = vapply(rnorm(n=sc,b,sqrt(b/s)),FUN=function(x) max(x,0.000001),1)

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

rgl2 = rg(1/9); rg13 = rg(1/3);      rgm = rg(1);     #rgn = rg(1/1000);
rLNl2 = rLN(1/9); rLN13 = rLN(1/3);  rLNm = rLN(1);  #rlNn = rLN(1/1000);
rCNl2 = rCN(1/9); rCN13 = rCN(1/3);  rCNm = rCN(1); #rCNn = rCN(1/1000);
rex = rexp();

d=rbind(
  #data.frame(m="LSD",r="Strict",n=lsd(be,b,1),v=0),
  #data.frame(m="LogDate",r="Strict",n=logD(be,b,1),v=0),
  data.frame(m="LSD",r="Gamma",n=lsd(be,b,rgl2),v="Var: 1/9",o=be*rgl2),
  data.frame(m="LogDate",r="Gamma",n=logD(be,b,rgl2),v="Var: 1/9",o=be*rgl2),
  data.frame(m="LSD",r="LogNormal",n=lsd(be,b,rLNl2),v="Var: 1/9",o=be*rLNl2),
  data.frame(m="LogDate",r="LogNormal",n=logD(be,b,rLNl2),v="Var: 1/9",o=be*rLNl2),   
  data.frame(m="LSD",r="LogUnif",n=lsd(be,b,rCNl2),v="Var: 1/9",o=be*rCNl2),
  data.frame(m="LogDate",r="LogUnif",n=logD(be,b,rCNl2),v="Var: 1/9",o=be*rCNl2), 
  
  data.frame(m="LSD",r="Gamma",n=lsd(be,b,rg13),v="Var:1/3",o=be*rg13),
  data.frame(m="LogDate",r="Gamma",n=logD(be,b,rg13),v="Var:1/3",o=be*rg13),
  data.frame(m="LSD",r="LogNormal",n=lsd(be,b,rLN13),v="Var:1/3",o=be*rLN13),
  data.frame(m="LogDate",r="LogNormal",n=logD(be,b,rLN13),v="Var:1/3",o=be*rLN13),
  data.frame(m="LSD",r="LogUnif",n=lsd(be,b,rCN13),v="Var:1/3",o=be*rCN13),
  data.frame(m="LogDate",r="LogUnif",n=logD(be,b,rCN13),v="Var:1/3",o=be*rCN13),
  
  
  data.frame(m="LSD",r="Gamma",n=lsd(be,b,rgm),v="Var: 1",o=be*rgm),
  data.frame(m="LogDate",r="Gamma",n=logD(be,b,rgm),v="Var: 1",o=be*rgm),
  data.frame(m="LSD",r="LogNormal",n=lsd(be,b,rLNm),v="Var: 1",o=be*rLNm),
  data.frame(m="LogDate",r="LogNormal",n=logD(be,b,rLNm),v="Var: 1",o=be*rLNm),
  data.frame(m="LSD",r="LogUnif",n=lsd(be,b,rCNm),v="Var: 1",o=be*rCNm),
  data.frame(m="LogDate",r="LogUnif",n=logD(be,b,rCNm),v="Var: 1",o=be*rCNm),
  
  data.frame(m="LSD",r="Exponential",n=lsd(be,b,rex),v="Var: 2.1",o=be*rex),
  data.frame(m="LogDate",r="Exponential",n=logD(be,b,rex),v="Var: 2.1",o=be*rex),
  
  data.frame(m="LSD",r="Strict",n=lsd(be,b,1),v="Var: 0",o=be*1),
  data.frame(m="LogDate",r="Strict",n=logD(be,b,1),v="Var: 0",o=be*1) 
);

d$n2=d$n^2

ggplot(aes(x=n), data=d)+
      geom_density(aes(color=m))+geom_histogram(aes(fill=m,y=..density..),position = "identity",binwidth = 0.05,alpha=0.4)+
  theme_classic()+facet_wrap(~interaction(r,v,sep=": "),scales="free",ncol=3)+ 
  scale_color_brewer(name="",palette = "Dark2")+scale_fill_brewer(name="",palette = "Dark2")+xlab("Penalty (without square)")+xlim(-10,10)+
  geom_vline(xintercept = 0,linetype=3)+theme(legend.position = c(0.87,0.2))+coord_cartesian(xlim=c(-3.3,3.3))
ggsave("product.pdf",width=7,height = 9)


# ggplot(aes(x=n^2), data=d[d$m=="LSD",])+
#   #geom_density(aes(color=m),kernel="b",adjust=2)+
#   geom_histogram(aes(fill=m),position = "identity",binwidth = 0.01,alpha=0.4)+
#   theme_classic()+facet_wrap(v~r,scales="free")+ 
#   scale_color_brewer(name="",palette = "Dark2")+scale_fill_brewer(name="",palette = "Dark2")+
#   xlab("Penalty")+xlim(0,100)+coord_cartesian(xlim=c(0,1))+#scale_x_log10(lim=c(0.0000001,1000))
#   geom_vline(xintercept = 0,linetype=3)+theme(legend.position = c(0.87,0.2))
# ggsave("product-square.pdf",width=10,height = 6)


d$oc= cut(log(d$o), (-2500:200)/50)
a = recast(oc+m+r+v~.,data=d,measure.var ="n",fun.aggregate=length)
a$n2=recast(oc+m+r+v~.,data=d,measure.var ="n2",fun.aggregate=median)$.
a$n=recast(oc+m+r+v~.,data=d,measure.var ="n",fun.aggregate=median)$.
a$c = a$.

ggplot(aes(y=c/sc,x=n2,color=m), data=a)+
  #geom_density(aes(color=m),kernel="b",adjust=2)+
  geom_point(alpha=0.6)+
  theme_classic()+theme(legend.position = c(0.8,0.15))+
  facet_wrap(~interaction(r,v,sep=": "),scales="free",ncol=3)+ 
  scale_color_brewer(name="",palette = "Dark2")+scale_fill_brewer(name="",palette = "Dark2")+
  xlab("Penalty")+scale_x_sqrt(lim=c(0,10))+scale_y_continuous(name="Frequency of data",labels=percent)
ggsave("product-penalty.pdf",width=7,height = 9)

ggplot(aes(x=log(c/sc/5),y=n2,color=m), data=a[(a$c/sc/5)>0.001,])+
  #geom_density(aes(color=m),kernel="b",adjust=2)+
  geom_point(alpha=0.6,size=0.5)+
  theme_classic()+theme(legend.position = c(0.8,0.15))+
  facet_wrap(~interaction(r,v,sep=": "),scales="free",ncol=3)+ 
  scale_color_brewer(name="",palette = "Dark2")+scale_fill_brewer(name="",palette = "Dark2")+
  scale_y_continuous(name="Penalty")+scale_x_continuous(name="Log likelihood (empirical)")
ggsave("compound-ll.pdf",width=7,height = 9)

######################### Compound


t = 500000
per = 10
sc = t/per
s=500
ee  = function(bt) vapply(rnorm(n=t,bt,sqrt(bt/s)),FUN=function(x) max(x,0.00001),1) 
ee  = function(bt) {a=rnorm(n=t,bt,sqrt(bt/s));a[a>0]}

bg = function(v) { ee( rgamma(sc,1/v,1/v)*b)}
bexp = function() {ee (stats::rexp(sc,log(2))*b )}
bLN = function(v) {ee (rlnorm(sc,meanlog = 0,sdlog =sqrt(log(1/2 + 1/2 *sqrt(1 + 4 *v))))*b)}
bCN = function(v,a=0.1) {
  x=(1:60)/40;
  p=x[which.min(abs(c((exp(2*x)-exp(-2*x))/4/x-(exp(x)-exp(-x))^2/4/x^2)-v))];
  ee( exp(runif(sc,-p,p))*b)}

summary(bCN(1/3)); var(bCN(1/3));
summary(bLN(1/2)); var(bLN(1/2));
summary(bg(1/2)); var(bg(1/2));
summary(bexp()); var(bexp())

bgl2 = bg(1/9);   bg13 = bg(1/3);    bgm = bg(1);     #bgn = bg(1/1000);
bLNl2 = bLN(1/9); bLN13 = bLN(1/3);  bLNm = bLN(1);  #blNn = bLN(1/1000);
bCNl2 = bCN(1/9); bCN13 = bCN(1/3);  bCNm = bCN(1); #bCNn = bCN(1/1000);
bex = bexp();
bs=ee(1*b)

lsd = function(be,b) {b/be-1}
logD = function(be,b) {log(b/be)}

d=rbind(
  #data.frame(m="LSD",r="Strict",n=lsd(be,b,1),v=0),
  #data.frame(m="LogDate",r="Strict",n=logD(be,b,1),v=0),
  data.frame(m="LSD",r="Gamma",n=lsd(bgl2,b),v="Var: 1/9",o=bgl2),
  data.frame(m="LogDate",r="Gamma",n=logD(bgl2,b),v="Var: 1/9",o=bgl2),
  data.frame(m="LSD",r="LogNormal",n=lsd(bLNl2,b),v="Var: 1/9",o=bLNl2),
  data.frame(m="LogDate",r="LogNormal",n=logD(bLNl2,b),v="Var: 1/9",o=bLNl2),   
  data.frame(m="LSD",r="LogUnif",n=lsd(bCNl2,b),v="Var: 1/9",o=bCNl2),
  data.frame(m="LogDate",r="LogUnif",n=logD(bCNl2,b),v="Var: 1/9",o=bCNl2), 
  
  data.frame(m="LSD",r="Gamma",n=lsd(bg13,b),v="Var:1/3",o=bg13),
  data.frame(m="LogDate",r="Gamma",n=logD(bg13,b),v="Var:1/3",o=bg13),
  data.frame(m="LSD",r="LogNormal",n=lsd(bLN13,b),v="Var:1/3",o=bLN13),
  data.frame(m="LogDate",r="LogNormal",n=logD(bLN13,b),v="Var:1/3",o=bLN13),
  data.frame(m="LSD",r="LogUnif",n=lsd(bCN13,b),v="Var:1/3",o=bCN13),
  data.frame(m="LogDate",r="LogUnif",n=logD(bCN13,b),v="Var:1/3",o=bCN13),
  
  
  data.frame(m="LSD",r="Gamma",n=lsd(bgm,b),v="Var: 1",o=bgm),
  data.frame(m="LogDate",r="Gamma",n=logD(bgm,b),v="Var: 1",o=bgm),
  data.frame(m="LSD",r="LogNormal",n=lsd(bLNm,b),v="Var: 1",o=bLNm),
  data.frame(m="LogDate",r="LogNormal",n=logD(bLNm,b),v="Var: 1",o=bLNm),
  data.frame(m="LSD",r="LogUnif",n=lsd(bCNm,b),v="Var: 1",o=bCNm),
  data.frame(m="LogDate",r="LogUnif",n=logD(bCNm,b),v="Var: 1",o=bCNm),
  
  data.frame(m="LSD",r="Exponential",n=lsd(bex,b),v="Var: 2.1",o=bex),
  data.frame(m="LogDate",r="Exponential",n=logD(bex,b),v="Var: 2.1",o=bex),
  
  data.frame(m="LSD",r="Strict",n=lsd(bs,b),v="Var: 0",o=bs),
  data.frame(m="LogDate",r="Strict",n=logD(bs,b),v="Var: 0",o=bs) 
);

ggplot(aes(x=n), data=d)+
  geom_density(aes(color=m))+geom_histogram(aes(fill=m,y=..density..),position = "identity",binwidth = 0.05,alpha=0.4)+
  theme_classic()+facet_wrap(~interaction(r,v,sep=": "),scales="free",ncol=3)+ 
  scale_color_brewer(name="",palette = "Dark2")+scale_fill_brewer(name="",palette = "Dark2")+xlab("Penalty (without square)")+xlim(-10,10)+
  geom_vline(xintercept = 0,linetype=3)+theme(legend.position = c(0.87,0.1))+coord_cartesian(xlim=c(-3.3,3.3))
ggsave("compound.pdf",width=7,height = 9)

#d$oc= cut(log(d$o), (-1000:50)/50)
#d$oc= cut(d$o, 50000)
#d$oc= cut(d$n,(-1000:1000000)/50)
dc = dcast(r+v+o~m,data=d,value.var = "n");
dc$oc  = cut(dc$LogDate,1000);
d = melt(dc,id.vars = c(1:3,6),value.name = "n",variable.name = "m");
d$n2=d$n^2;

a = recast(oc+m+r+v~.,data=d,measure.var ="n",fun.aggregate=length)
a$n2=recast(oc+m+r+v~.,data=d,measure.var ="n2",fun.aggregate=median)$.
a$n=recast(oc+m+r+v~.,data=d,measure.var ="n",fun.aggregate=median)$.
a$c = a$.

ggplot(aes(y=(c/t),x=n2,color=m), data=a[log(a$c/t)>-7&a$n2<1000,])+
  #geom_density(aes(color=m),kernel="b",adjust=2)+
  geom_point(alpha=0.6,size=0.5)+
  theme_classic()+theme(legend.position = c(0.8,0.15))+
  facet_wrap(~interaction(r,v,sep=": "),scales="free",ncol=3)+ 
  scale_color_brewer(name="",palette = "Dark2")+scale_fill_brewer(name="",palette = "Dark2")+
  scale_x_continuous(name="Penalty")+scale_y_log10(name="Likelihood (empirical)")
ggsave("compound-penalty.pdf",width=7,height = 9)

ggplot(aes(x=log(c/t),y=n2,color=m), data=a[log(a$c/t)>-7&a$n2<1000,])+
  #geom_density(aes(color=m),kernel="b",adjust=2)+
  #geom_line(aes(group=interaction(m,n>0)),size=0.4)+
  geom_point(alpha=0.5,size=0.5)+
  theme_classic()+theme(legend.position = c(0.8,0.15))+
  facet_wrap(~interaction(r,v,sep=": "),scales="free",ncol=3)+ 
  scale_color_brewer(name="",palette = "Dark2")+scale_fill_brewer(name="",palette = "Dark2")+
  scale_y_log10(name="Penalty")+scale_x_continuous(name="Log likelihood (empirical)")
ggsave("compound-ll.pdf",width=7,height = 9)
