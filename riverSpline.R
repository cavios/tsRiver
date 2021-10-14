require(TMB)
source('tsRiverLib.R')

## ======================= Parameterization of phi
f <- function(x) 2/(1 + exp(-2 * x)) - 1
invf <- function(y) -0.5 * log(2/(y + 1) - 1)




## ======================= Fit model

compile("riverSpline.cpp")
dyn.load(dynlib("riverSpline"))

dat<-read.table('tsRiver_missouri_csssi.dat',header=TRUE)
# here we select a 200km reach
id<-which(dat$dist > VS-100000 & VS+100000)
dat<-dat[id,]

VS<-220672.7
VSname<-'sioux'



o<-order(dat$time)
dato<-dat[o,]
startH<-mean(dat$height[order(dat$dist)[1:100]])


#define values for time steps and number of nodes in the spline functions
tdim <- 500
xdim <- 7
tidx <- assignPretty(dato$time, tdim)
myrange<-range(dato$dist);xknot<-seq(myrange[1],myrange[2],length=xdim)

#define the number of iterations
nit<-5

# data transfered to the cpp part
H<-dato$height
data <- list(H=H,
             x=dato$dist,
             w=0*H+1,
             xknot=xknot,
             VS=VS,
             tidx=tidx$ct,
             satid=dato$satId
             )

# initial parameter values tranfered to cpp part
parameters<- list(eta=rep(0,tdim),
                      mu=c(startH,rep(0,length(data$xknot)-1)),
                      logSd=rep(0,length(unique(data$satid))),
                      logAscale=rep(0,length(data$xknot)),
                      bias=rep(0,length(unique(data$satid))-1),
                      transf_phi=rep(invf(0.5),1)
                      )

allnll<-numeric(nit)
#first iteration
obj <- MakeADFun(data,parameters,random="eta",DLL="riverSpline")
#minimization
opt<-nlminb(obj$par,obj$fn,obj$gr,control=list(eval.max=2000, iter.max=2000))
allnll[1]<-opt$objective
mypred<-obj$report()$pred
mypredsd<-obj$report()$predsd


# recompute weights
for(i in 2:length(allnll)){
    cutp <- 0.01
    tailp <- 1-pnorm(abs(data$H-mypred), mean=0, sd=mypredsd)
    w <- ifelse(tailp<cutp, tailp/cutp, 1)
    eps <- 1.0e-16
    data$w= ifelse(w<eps, eps, w)
    obj <- MakeADFun(data,parameters,random="eta",DLL="riverSpline")
    #minimization
    opt<-nlminb(obj$par,obj$fn,obj$gr, control=list(eval.max=3000, iter.max=3000))
    allnll[i]<-opt$objective
    mypred<-obj$report()$pred
    mypredsd<-obj$report()$predsd
}

rep<-sdreport(obj)
myH<-as.list(rep,"Est",report=TRUE)$hVS
mySD<-as.list(rep,"Std",report=TRUE)$hVS
time<-tidx$xx
out<-data.frame(time=time,H=myH,Hsd=mySD)
ofile<-paste0('riverTS_',VSname,'_',formatC(tdim,width=4,flag="0"),'_',formatC(xdim,width=2,flag="0"),'_',formatC(nit,width=2,flag="0"),'.dat',sep='')
write.table(out,ofile, row.names=FALSE, quote=FALSE)



#plot of time series 
insitu<-read.table('siouxD.txt',header=TRUE)
xlim1<-c(2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020,2021)
xlim2<-c(2011,2012,2013,2014,2015,2016,2017,2018,2019,2020,2021,2022)

par(mfrow=c(2,6))

for(i in 1:length(xlim1)){
plot(out[,1],out[,2],t='l',xlim=c(xlim1[i],xlim2[i]),col='blue',lwd=4)
lines(insitu$time,insitu$elev,col='red',lwd=4)
lines(out[,1],out[,2]-2*out[,3],col='blue',lty=2)
lines(out[,1],out[,2]+2*out[,3],col='blue',lty=2)

}

   

