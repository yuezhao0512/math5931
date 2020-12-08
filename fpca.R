data_1=read.csv(file.choose(),header = FALSE)
data=data_1[,-1]
data = data.matrix(data)
data = t(data)

library(fda)
daybasis365 = create.fourier.basis(c(0,246),246)
harmLfd = vec2Lfd(c(0,(2*pi/246)^2,0), c(0, 246))
earnfdPar = fdPar(daybasis365,harmLfd,1e4)
earnfd = smooth.basis(1:246,data,earnfdPar)

#quartz()
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
plot(earnfd$fd,xlab='day',ylab='Earning Ratio',cex.lab=1.5,cex.axis=1.5)

#daily$place
#length(daily$place)

earnvar = var.fd(earnfd$fd)
tvvals = eval.bifd(1:246,1:246,earnvar)

par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
contour(1:246,1:246,tvvals,xlab='day',ylab='day',cex.lab=1.5,cex.axis=1.5)
#install.packages("fields")
library(fields)

image.plot(1:246,1:246,tvvals,xlab='day',ylab='day',cex.lab=1.5,cex.axis=1.5)
earn.cor = cor.fd(1:246,earnfd$fd)
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
contour(1:246,1:246,earn.cor,xlab='day',ylab='day',cex.lab=1.5,cex.axis=1.5)
image.plot(1:246,1:246,earn.cor,xlab='day',ylab='day',cex.lab=1.5,cex.axis=1.5)




earnpca = pca.fd(earnfd$fd,nharm=20)
names(earnpca)
earnpca$varprop
#temppca$values are the eigenvalues
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
plot(earnpca$values[1:20],xlab='component',ylab='variance',col="red",
     cex.lab=1.5,cex.axis=1.5,cex=2)

par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
plot(cumsum(earnpca$values[1:20])/sum(earnpca$values),xlab='Number of Components',
     ylab='cumulative variance explained',col=2,cex.lab=2,
     cex.axis=2,cex=2)
abline(h=0.90)


harmfd = earnpca$harmonics
harmvals = eval.fd(1:246,harmfd)
dim(harmvals) # The top 20 FPCs


par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
plot(1:246,harmvals[,1],xlab='day',ylab='PC1',
     lwd=4,lty=1,cex.lab=2,cex.axis=2,type='l')
title('First Principle Component Functions')
abline(h=0, col='red')


par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
plot(1:246,harmvals[,2],xlab='day',ylab='PCs2',
     lwd=4,lty=1,cex.lab=2,cex.axis=2,type='l')
title('Second Principle Component Functions')
abline(h=0, col='red')


par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
plot(1:246,harmvals[,3],xlab='day',ylab='PCs3',
     lwd=4,lty=1,cex.lab=2,cex.axis=2,type='l')
abline(h=0, col='red')

par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
plot(1:246,harmvals[,4],xlab='day',ylab='PCs4',
     lwd=4,lty=1,cex.lab=2,cex.axis=2,type='l')
abline(h=0, col='red')

par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
plot(1:246,harmvals[,5],xlab='day',ylab='PCs5',
     lwd=4,lty=1,cex.lab=2,cex.axis=2,type='l')
abline(h=0, col='red')

par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
plot(1:246,harmvals[,6],xlab='day',ylab='PCs6',
     lwd=4,lty=1,cex.lab=2,cex.axis=2,type='l')
abline(h=0, col='red')

par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
plot(1:246,harmvals[,7],xlab='day',ylab='PCs7',
     lwd=4,lty=1,cex.lab=2,cex.axis=2,type='l')
abline(h=0, col='red')

par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
plot(1:246,harmvals[,8],xlab='day',ylab='PCs8',
     lwd=4,lty=1,cex.lab=2,cex.axis=2,type='l')
abline(h=0, col='red')

par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
plot(1:246,harmvals[,9],xlab='day',ylab='PCs9',
     lwd=4,lty=1,cex.lab=2,cex.axis=2,type='l')
abline(h=0, col='red')

par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
matplot(1:246,harmvals[,1:6],xlab='day',ylab='PCs',
        lwd=4,lty=1,cex.lab=2.5,cex.axis=2.5,type='l')
legend(0,-0.04,c('PC1','PC2','PC3','PC4','PC5','PC6'),col=1:4,lty=1,lwd=5)
title('Functional Principle Component Functions')

par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
plot(temppca$scores[,1:2],xlab='PC Score 1',ylab='PC Score 2',col=4,
     cex.lab=1.5,cex.axis=1.5,cex=1)
text(earnpca$scores[,1],earnpca$scores[,2],labels=data_1[,1],cex=1)
title('First FPC score VS Second FPC score')

par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
plot(earnpca$scores[,2:3],xlab='PC Score 2',ylab='PC Score 3',col=4,
     cex.lab=1.5,cex.axis=1.5,cex=1)
text(earnpca$scores[,2],earnpca$scores[,3],labels=data_1[,1],cex=1)


par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
plot(earnpca$scores[,3:4],xlab='PC Score 3',ylab='PC Score 4',col=4,
     cex.lab=1.5,cex.axis=1.5,cex=1)
text(earnpca$scores[,3],earnpca$scores[,4],labels=data_1[,1],cex=1)







