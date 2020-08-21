
library(MASS)
library(sampling)

calib.logist=function (Xs, d, total,max_iter = 500) 
{
    q = rep(1, length(d))
    EPS = .Machine$double.eps
    EPS1 = 1e-06
    n = length(d)  
           
       lambda = as.matrix(rep(0, ncol(Xs)))
        w1 = as.vector(d * (1+exp((Xs%*%lambda*q))))
        for (l in 1:max_iter) {
            phi = t(Xs) %*% w1 - total
            T1 = t(Xs * w1)
            phiprim = T1 %*% Xs
            lambda = lambda - ginv(phiprim, tol = EPS) %*% phi
            w1 = as.vector(d * (1+exp((Xs %*% lambda * q))))
            if (any(is.na(w1)) | any(is.infinite(w1))) {
                warning("No convergence")
                g = NA
                break
            }
            tr = crossprod(Xs, w1)
            if (max(abs(tr - total)/total) < EPS1) 
                break
        }
        if (l == max_iter) {
            warning("No convergence")
            g = NA
        }
        else g = w1/d
  
        g
}
#######################
calib.log.c=function (Xs, d, total,bounds=c(l,u,c), max_iter = 500) 
{
   q = rep(1, length(d))
    EPS = .Machine$double.eps
    EPS1 = 1e-06
    n = length(d)
    lambda = as.matrix(rep(0, n))
    lambda1 = ginv(t(Xs * d * q) %*% Xs, tol = EPS) %*% (total - 
        as.vector(t(d) %*% Xs))
  
        A = (bounds[2] - bounds[1])/((bounds[3] - bounds[1]) * (bounds[2] - bounds[3]))
        u = rep(1, length(d))
        F = (bounds[1]*(bounds[2] - bounds[3])+bounds[2]*(bounds[3]-bounds[1]) *u)/(bounds[2]-bounds[3]+(bounds[3]-bounds[1])*u)
        w1 = as.vector(d * F)
        T = t(Xs * w1)
        phiprim = ginv(T %*% Xs, tol = EPS)
        g = F
        tr = crossprod(Xs, w1)
        if (max(abs(tr - total)/total) > EPS1 | any(g < bounds[1]) | 
            any(g > bounds[2])) {
            lambda1 = rep(0, ncol(Xs))
            list = 1:length(g)
            t2 = total
            Xs1 = Xs
            d1 = d
            g1 = g
            q1 = q
            list1 = 1:length(g)
            for (l in 1:max_iter) {
                if (any(g < bounds[1]) | any(g > bounds[2])) {
                  g[g < bounds[1]] = bounds[1]
                  g[g > bounds[2]] = bounds[2]
                  list = (1:length(g))[g > bounds[1] & g < bounds[2]]
                  if (length(list) != 0) {
                    g1 = g[list]
                    t2 = total - as.vector(t(g[-list] * d[-list]) %*% 
                      Xs[-list, ])
                    Xs1 = Xs[list, ]
                    d1 = d[list]
                    q1 = q[list]
                    list1 = list
                  }
                }
                t1 = as.vector(t(d1) %*% Xs1)
                phi = t(Xs1) %*% as.vector(d1 * g1) - t1
                T = t(Xs1 * as.vector(d1 * g1))
                phiprime = T %*% Xs1
                lambda1 = lambda1 - ginv(phiprime, tol = EPS) %*% 
                  (as.vector(phi) - t2 + t1)
                u = exp(A * (Xs1 %*% lambda1 * q1))
                F = g1 = (bounds[1] * (bounds[2] - bounds[3]) + bounds[2] * 
                  (bounds[3] - bounds[1]) * u)/(bounds[2] - bounds[3] + (bounds[3]- 
                  bounds[1]) * u)
                if (any(is.na(g1))) {
                  warning("no convergence")
                  g1 = g = NA
                  break
                }
                g[list1] = g1
                tr = crossprod(Xs, g * d)
                if (max(abs(tr - total)/total) < EPS1 & all(g >= 
                  bounds[1] & g <= bounds[2])) 
                  break
            }
            if (l == max_iter) {
                cat("no convergence in", max_iter, "iterations with the given bounds. \n")
                cat("the bounds for the g-weights are:", min(g), 
                  " and ", max(g), "\n")
                cat(" and the g-weights are given by g\n")
                g = NA
            }
        }
    g
}
####################################
N=10000
n=200
B=1000
teta=3

T.Lin.sam=T.Log.sam=T.Log.c.sam=T.Rak.sam=T.Lin.pop=T.Log.pop=T.Log.c.pop=T.Rak.pop=0
pb=winProgressBar(title = "progress bar", min = 0,max = B, width = 300)

for(i in 1:B){
z1=rexp(N)
e=rnorm(N)
y=2+z1+e
s=sample(1:N,n,replace=FALSE)
p=1/(1+0.2*exp(2*(log(z1[s])+0.5)))

r=numeric()
for(a in 1:n){r[a]=rbinom(1,1,p[a])}

y.S=y[s]
y.k=y.S[r==1]

pi.k=n/N
d.k=1/pi.k
####################################
#To the sample with xk=(1 log(zk))

x.U=matrix(c(rep(1,N),log(z1)),ncol=2)

x.S=x.U[s,]
x.R=x.S[r==1,]

sjx.U=colSums(x.U)
snx.R=d.k*colSums(x.R)
snx.S=d.k*colSums(x.S)

w.Lin.sam=calib(x.R,d=1/rep(pi.k,nrow(x.R)),snx.S,method="linear")
w.Log.sam=calib.logist(x.R,d=1/rep(pi.k,nrow(x.R)),snx.S,max_iter =500) 
w.Rak.sam=calib(x.R,d=1/rep(pi.k,nrow(x.R)),snx.S,method="raking")
w.Log.c.sam=calib.log.c(x.R,d=1/rep(pi.k,nrow(x.R)),snx.S,bounds=c(l=1,u=1000000,c=2),max_iter =500) 

T.Lin.sam[i]=sum(w.Lin.sam*y.k)/sum(w.Lin.sam)
T.Rak.sam[i]=sum(w.Rak.sam*y.k)/sum(w.Rak.sam)
T.Log.sam[i]=sum(w.Log.sam*y.k)/sum(w.Log.sam)
T.Log.c.sam[i]=sum(w.Log.c.sam*y.k)/sum(w.Log.c.sam)
############################
#To the population with xk=(1 zk)

x.U=matrix(c(rep(1,N),z1),ncol=2)

x.S=x.U[s,]
x.R=x.S[r==1,]

sjx.U=colSums(x.U)
snx.R=d.k*colSums(x.R)
snx.S=d.k*colSums(x.S)

w.Lin.pop=calib(x.R,d=1/rep(pi.k,nrow(x.R)),sjx.U,method="linear")
w.Log.pop=calib.logist(x.R,d=1/rep(pi.k,nrow(x.R)),sjx.U,max_iter =500) 
w.Rak.pop=calib(x.R,d=1/rep(pi.k,nrow(x.R)),sjx.U,method="raking")
w.Log.c.pop=calib.log.c(x.R,d=1/rep(pi.k,nrow(x.R)),sjx.U,bounds=c(l=1,u=10000000,c=2),max_iter =500) 

T.Lin.pop[i]=sum(w.Lin.pop*y.k)/sum(w.Lin.pop)
T.Rak.pop[i]=sum(w.Rak.pop*y.k)/sum(w.Rak.pop)
T.Log.pop[i]=sum(w.Log.pop*y.k)/sum(w.Log.pop)
T.Log.c.pop[i]=sum(w.Log.c.pop*y.k)/sum(w.Log.c.pop)

setWinProgressBar(pb,i, title=paste( round(i/B*100, 0),"%","B=",i))}
close(pb)
T.hat=matrix(c(T.Lin.sam,T.Rak.sam,T.Log.sam,T.Log.c.sam,T.Lin.pop,T.Rak.pop,T.Log.pop,T.Log.c.pop),ncol=8)
T.hat=na.omit(T.hat)

Bias=apply(T.hat,2,mean)-teta
Var=apply(T.hat,2,var)
R.Bias=Bias/teta
R.MSE=sqrt((Var+Bias^2))/(teta^2)

names(R.Bias)=names(R.MSE)=names(Var)=c("T.Lin.sam","T.Rak.sam","T.Log.sam","T.Log.c.sam","T.Lin.pop","T.Rak.pop","T.Log.pop","T.Log.c.pop")

R.Bias*100
R.MSE*10000
Var*1000


#################################################################
data_as=as.data.frame(read.table(file.choose(), skip=0, sep=",", header=T))
attach(data_as)
head(data_as)

N=nrow(data_as)
n=300
B=1000
dbp=0.10
teta=mean(data_as$Daramad/10000000)

T.Lin.sam=T.Log.sam=T.Rak.sam=T.Log.c=sam.w.lin=sam.w.rak=sam.w.log=sam.w.L.U=test.Lin=test.Log=test.Rak=
test.Log.L.U=NA
pb=winProgressBar(title="progress bar", min = 0,max = B, width = 300)

for(i in 1:B){

s=sample(1:N,n)
data1=data_as[s,]
o=order(data1$NHazineh)
data=data1[o,]
L1=dbp*n*.99
L2=dbp*n*.01
r=sample(c(rep(0,L1),rep(1,(n-L1-L2)),rep(0,L2)))

y.S=data$Daramad/10000000
y.k=y.S[r==1]

pi.k=n/N
d.k=1/pi.k

x.U=matrix(c(rep(1,N),data_as$C01,data_as$NHazineh/10000000),ncol=3)
x.S=matrix(c(rep(1,n),data$C01,data$NHazineh/10000000),ncol=3)
x.R=x.S[r==1,]

sjx.U=colSums(x.U)
snx.R=d.k*colSums(x.R)
snx.S=d.k*colSums(x.S)

g.sam=c(t(sjx.U-snx.R)%*%solve(d.k*t(x.R)%*%x.R))
w.Lin.sam=d.k*(1+rowSums(t(g.sam*t(x.R))))
#w.Lin.sam=calib(x.R,d=1/rep(pi.k,nrow(x.R)),sjx.U,method="linear")
w.Log.sam=calib.logist(x.R,d=1/rep(pi.k,nrow(x.R)),sjx.U,max_iter =500) 
w.Rak.sam=calib(x.R,d=1/rep(pi.k,nrow(x.R)),sjx.U,method="raking")
w.Log.c.sam=calib.log.c(x.R,d=1/rep(pi.k,nrow(x.R)),sjx.U,bounds=c(l=0,u=6,c=1),max_iter =100) 

T.Lin.sam[i]=sum(w.Lin.sam*y.k)/sum(w.Lin.sam)
T.Rak.sam[i]=sum(w.Rak.sam*y.k)/sum(w.Rak.sam)
T.Log.sam[i]=sum(w.Log.sam*y.k)/sum(w.Log.sam)
T.Log.c.sam[i]=sum(w.Log.c.sam*y.k)/sum(w.Log.c.sam)
sam.w.lin[i]=sum(w.Lin.sam)
sam.w.rak[i]=d.k*sum(w.Rak.sam)
sam.w.log[i]=sum(w.Log.sam)
sam.w.L.U[i]=d.k*sum(w.Log.c.sam)

test.Lin[i]=sum(w.Lin.sam*x.R)-sum(sjx.U)
test.Log[i]=sum(w.Log.sam*x.R)-sum(sjx.U)
test.Rak[i]=(d.k*sum(w.Rak.sam*x.R))-sum(sjx.U)
test.Log.L.U[i]=(d.k*sum(w.Log.c.sam*x.R))-sum(sjx.U)


setWinProgressBar(pb,i, title=paste( round(i/B*100, 0),"%","B=",i))}
close(pb)
T.hat=matrix(c(T.Lin.sam,T.Rak.sam,T.Log.sam,T.Log.c.sam),ncol=4)
T.hat=na.omit(T.hat)
sum=matrix(c(sam.w.lin,sam.w.rak,sam.w.log,sam.w.L.U),ncol=4)
test=matrix(c(test.Lin,test.Rak,test.Log,test.Log.L.U),ncol=4)

sum1=apply(sum,2,mean)-N
test1=apply(test,2,mean)

Bias=apply(T.hat,2,mean)-teta
Var=apply(T.hat,2,var)
R.Bias=Bias/teta
R.MSE=sqrt((Var+Bias^2))/teta
names(R.Bias)=names(R.MSE)=names(Var)=c("T.Lin","T.Rak","T.Log","T.Log.c")
names(sum1)=c("Lin","Rak","Log","Log.L.U")
names(test1)=c("Lin","Rak","Log","Log.L.U")

R.Bias
R.MSE
sum1
test1
#-------------------------
n=c(100,200,500)
logit=c(0.09380140,0.07147505,0.06269165)
raking=c(0.09402484,0.07163095,0.06272275)
GREG=c(0.09336542,0.07122446,0.06261960)
ylim=range(logit,raking,GREG)

plot(n,logit,type="c",pch=19,lty=2,xlab="sample size",ylab="ralative efficiently",xaxt = "n",ylim=ylim)
axis(1, n)
points(n,raking,type="o",pch=15,lty=3)
points(n,GREG,type="b",pch=12,lty=4)
legend("bottomright",c("Logit","Raking","GREG"),bty="sample size",lty=c(2,3,4),pch =c(19,15,12))
