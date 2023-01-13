
library(CVAS)
 df=data.frame(product=rep(c("a","b","c"),each=10),
 subject=rep(paste0("S",1:10),3),acide=rnorm(30),amer=rnorm(30))
#  respca=PCA(df)
#  calculateEllipses(df,respca$EigenVectors)
#  # test verification:
#  centre=c(coorCentre[i,1],coorCentre[i,2])
#  # n*t(g-mu)(VarTheorique)^(-1)(g-mu) suit une loi du chi 2 (Saporta)
# # # et donc on a t(g-mu)(VarTheorique)^(-1)(g-mu) < chi2(p , 1-alpha)/n est l'equation de l'ellipsoide de confiance pour mu.
# # # on suppose ici que VarTheorique=matCov
# # # on a alors l'ellipse suivante qui correspond ? ce qu'on obtient en mettant
# # covarianceMatrix2=matCov
# # # c=-2*qf(confInt,2,nb.suj-2)/(nb.suj-2)
#
#
# # # en utilisant la matrice de covariance divisee par n
# # covarianceMatrix2=matCov2
# # c=-(qchisq(confInt,2))
#
# # # c=-2*qf(confInt,2,nb.suj-2)/(nb.suj-2)
#
#  saportaEllipses=function(covarianceMatrix2,centre,nb.suj,test="f",confInt=0.95)
#  {
#  if(test=="chi"){c=-(qchisq(confInt,2))/nb.suj;print(-c)}
#  if(test=="f"){ c=-2*qf(confInt,2,nb.suj-2)/(nb.suj-2) ;print(-c)}
#  if(test=="sas"){ c=-(2/(nb.suj-2))*qf(confInt,2,nb.suj-2)*(nb.suj-1)/nb.suj ;print(-c)}
#  vap=eigen(solve(covarianceMatrix2))
#  vep2=eigen(solve(covarianceMatrix2))$vectors
#  a=sqrt(-c/(vap$values[1]))
#  b=sqrt(-c/(vap$values[2]))
#  x=rep(NA,length(seq(0,2*pi,0.1)));y=rep(NA,length(seq(0,2*pi,0.1)))
#  for (j in 1:length(seq(0,2*pi,0.1)))
#  {
#  x[j]= a*cos(seq(0,2*pi,0.1)[j])
#  y[j]= b*sin(seq(0,2*pi,0.1)[j])
#  }
#  covariancemat=vep2%*%diag(vap$values)%*%t(vep2)
#  x2=cbind(x,y)%*%t(as.matrix(vep2))[,1]
#  y2=cbind(x,y)%*%t(as.matrix(vep2))[,2]
#  x3=centre[1]+x2
#  y3=centre[2]+y2
#  l=cbind(x3,y3)
#  return(l)
#  }
#
#
#  l_chi=saportaEllipses(covarianceMatrix=matCov,nb.suj=nb.suj,test="chi")
#  l_f=saportaEllipses(matCov,nb.suj=nb.suj,test="f")
#  l_sas=saportaEllipses(matCov,nb.suj=nb.suj,test="sas")
#
# plot(ellipse::ellipse(x=covarianceMatrix,centre=c(coorCentre[i,1],coorCentre[i,2]),level=confInt,t=quant2),type="l",col="pink")
#
#  lines(x=l_f[,1],y=l_f[,2],type="l",col="red")
#  lines(x=l_chi[,1],y=l_chi[,2],col="blue",type='l')
#  lines(x=l_sas[,1],y=l_sas[,2],col="orange",type='l')
# lines(ellipse::ellipse(x=covarianceMatrix,centre=c(coorCentre[i,1],coorCentre[i,2]),level=confInt,t=quant),type="l",col="green")
#
# lines(ellipse::ellipse(x=covarianceMatrix,centre=c(coorCentre[i,1],coorCentre[i,2]),level=confInt,t=quant_sas),type="l",col="yellow")
#
#
# product=rep(c("a","b","c"),each=10);
# subject=rep(paste0("S",1:10),3);acide=rnorm(30);amer=rnorm(30)
# X=matrix(NA,30,2);X[,1]=acide;X[,2]=amer
#  df=data.frame(product=product,subject=subject,acide=acide,amer=amer)
#  resplsda=plsda(X=X,Y=product)
#  res=plotIndiv(resplsda,ellipse=TRUE)
# # respca=PCA(df)
#  resEll=calculateEllipses(df,resplsda$loadings,ellipseType="barycentric",confInt=0.9)
#  resEll2=calculateEllipses(df,resplsda$loadings,ellipseType="individual",confInt=0.9)
#  resEll3=calculateEllipses(df,resplsda$loadings,ellipseType="individual",confInt=0.9)
#
#  plot(res$df.ellipse[,1:2],xlim=c(min(res$df.ellipse[,c(1,3,5)]),max(res$df.ellipse[,c(1,3,5)])),
#         ylim=c(min(res$df.ellipse[,c(2,4,6)]),max(res$df.ellipse[,c(2,4,6)])),type="l",col="blue",lwd=2)
#  lines(res$df.ellipse[,3:4],type="l",col="orange",lwd=2)
#  lines(res$df.ellipse[,5:6],type="l",col="grey",lwd=2)
#  col=c("blue","orange","grey")
#  points(resEll$indivCoord$V1,resEll$indivCoord$V2,col=col[as.factor(resEll$indivCoord$product)],pch=16)
#  points(resEll$centers$V1,resEll$centers$V2,col=col[as.factor(resEll$centers$product)],pch=16,cex=2)
# lines(resEll$ellPoints[resEll$ellPoints[,"product"]=="a","x"],resEll$ellPoints[resEll$ellPoints[,"product"]=="a","y"],col=col[1])
# lines(resEll$ellPoints[resEll$ellPoints[,"product"]=="b","x"],resEll$ellPoints[resEll$ellPoints[,"product"]=="b","y"],col=col[2])
# lines(resEll$ellPoints[resEll$ellPoints[,"product"]=="c","x"],resEll$ellPoints[resEll$ellPoints[,"product"]=="c","y"],col=col[3])
# lines(resEll2$ellPoints[resEll$ellPoints[,"product"]=="a","x"],resEll2$ellPoints[resEll$ellPoints[,"product"]=="a","y"],col=col[1])
# lines(resEll2$ellPoints[resEll$ellPoints[,"product"]=="b","x"],resEll2$ellPoints[resEll$ellPoints[,"product"]=="b","y"],col=col[2])
# lines(resEll2$ellPoints[resEll$ellPoints[,"product"]=="c","x"],resEll2$ellPoints[resEll$ellPoints[,"product"]=="c","y"],col=col[3])
# lines(resEll3$ellPoints[resEll$ellPoints[,"product"]=="a","x"],resEll3$ellPoints[resEll$ellPoints[,"product"]=="a","y"],col=col[1],lty=2)
# lines(resEll3$ellPoints[resEll$ellPoints[,"product"]=="b","x"],resEll3$ellPoints[resEll$ellPoints[,"product"]=="b","y"],col=col[2],lty=2)
# lines(resEll3$ellPoints[resEll$ellPoints[,"product"]=="c","x"],resEll3$ellPoints[resEll$ellPoints[,"product"]=="c","y"],col=col[3],lty=2)
#
#
# plot(seq(0, 1, length.out=100),(1/n)*qchisq(seq(0, 1, length.out=100),df=2),type="l")
# lines(seq(0, 1, length.out=100),(2*(n-1)/(n-2))*qf(seq(0, 1, length.out=100),df1=2,df2=n-2),col="red")
# lines(seq(0, 1, length.out=100),2*qchisq(seq(0, 1, length.out=100),df=2),col="green")
#
# n=3
# plot(seq(0, 0.2, length.out=100),(1/n)*qf(seq(0, 1, length.out=100),df1=2,df2=2*n-3),type="l")
# lines(seq(0, 0.2, length.out=100),(2*n-3)/((n-1)*(n-2))*qf(seq(0, 1, length.out=100),df1=2,df2=n-2),col="red")
#
# #lines(seq(0, 1, length.out=100),2*qchisq(seq(0, 1, length.out=100),df=2),col="green")
