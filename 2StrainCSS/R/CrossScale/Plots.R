r20rates <- numeric(0)
for (i in 1:length(r20births))
{
  r20rates <- c(r20rates, 2.0)
}


plot(r20$V2, r20$V1, main="Particles Transmitted over Time, r=2.0", 
     xlab="Time", ylab="Particles Transmitted", pch=1, cex=0.5)

plot(r19$V2, r19$V1, main="Particles Transmitted over Time, r=1.9", 
     xlab="Time", ylab="Particles Transmitted", pch=1, cex=0.5)

plot(r18$V2, r18$V1, main="Particles Transmitted over Time, r=1.8", 
     xlab="Time", ylab="Particles Transmitted", pch=1, cex=0.5)

plot(r17$V2, r17$V1, main="Particles Transmitted over Time, r=1.7", 
     xlab="Time", ylab="Particles Transmitted", pch=1, cex=0.5)

plot(r16$V2, r16$V1, main="Particles Transmitted over Time, r=1.6", 
     xlab="Time", ylab="Particles Transmitted", pch=1, cex=0.5)

plot(r15$V2, r15$V1, main="Particles Transmitted over Time, r=1.5", 
     xlab="Time", ylab="Particles Transmitted", pch=1, cex=0.5)

plot(r14$V2, r14$V1, main="Particles Transmitted over Time, r=1.4", 
     xlab="Time", ylab="Particles Transmitted", pch=1, cex=0.5)

plot(r13$V2, r13$V1, main="Particles Transmitted over Time, r=1.3", 
     xlab="Time", ylab="Particles Transmitted", pch=1, cex=0.5)

plot(r12$V2, r12$V1, main="Particles Transmitted over Time, r=1.2", 
     xlab="Time", ylab="Particles Transmitted", pch=1, cex=0.5)

plot(r11$V2, r11$V1, main="Particles Transmitted over Time, r=1.1", 
     xlab="Time", ylab="Particles Transmitted", pch=1, cex=0.5)

plot(r10$V2, r10$V1, main="Particles Transmitted over Time, r=1.0", 
     xlab="Time", ylab="Particles Transmitted", pch=1, cex=0.5)

plot(r09$V2, r09$V1, main="Particles Transmitted over Time, r=0.9", 
     xlab="Time", ylab="Particles Transmitted", pch=1, cex=0.5)

plot(r08$V2, r08$V1, main="Particles Transmitted over Time, r=0.8", 
     xlab="Time", ylab="Particles Transmitted", pch=1, cex=0.5)

plot(r07$V2, r07$V1, main="Particles Transmitted over Time, r=0.7", 
     xlab="Time", ylab="Particles Transmitted", pch=1, cex=0.5)

plot(r06$V2, r06$V1, main="Particles Transmitted over Time, r=0.6", 
     xlab="Time", ylab="Particles Transmitted", pch=1, cex=0.5)

plot(r05$V2, r05$V1, main="Particles Transmitted over Time, r=0.5", 
     xlab="Time", ylab="Particles Transmitted", pch=1, cex=0.5)

`r05` <- read.csv("~/Documents/Courses/SIRModels/2StrainCSS/2StrainCSS/Build/Products/Debug/formatted.txt", header=F)

