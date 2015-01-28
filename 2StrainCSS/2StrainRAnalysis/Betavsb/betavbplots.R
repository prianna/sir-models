b <- c(b, 10^(-5.0))
beta <- c(beta, 2.28/T)

b10 <- read.csv("~/Documents/Courses/SIRModels/2StrainCSS/2StrainCSS/Build/Products/Debug/formatted.txt", header=F)
children <- numeric(0)
time <- numeric(0)


b10 <- read.csv("~/Documents/Courses/SIRModels/2StrainCSS/2StrainCSS/Build/Products/Debug/formatted.txt", header=F)
contact <- 0
nocontact <- 0
kids <- 0
for (i in 1:length(b10$V7))
{
  if ( b10$V7[i] > 0 && b10$V6[i] == 0)
    {
      nocontact = nocontact+1
    }
  if (b10$V7[i] > 0 && b10$V6[i] > 0)
    {
      contact = contact+1
    }
  if(b10$V5[i] == 0)
  {
    kids <- kids+b10$V6[i]
  }
}

prob10 <- nocontact/contact
R0 <- c(kids/10000, R0)


t4r20 <- time
b4r20 <- children

plot(time,children, xlim=x, main="Particles Transmitted over Time, b=1e-4", 
     xlab="Time", ylab="Particles Transmitted (Beta = 0.14)", pch=1, cex=0.5)

points(t1, b1, type="p", col="red")
points(t1r10, b1r10, type="p", col="blue")
points(t1r12, b1r12, type="p", col="green")
points(t1r14, b1r14, type="p", col="purple")
points(t1r18, b1r18, type="p", col="brown")

plot(t1r20, b1r20, xlim=x, main="Particles Transmitted over Time, b=1e-4", 
     xlab="Time", ylab="Particles Transmitted (Beta = 0.14)", pch=1, cex=0.5)

plot(t2r20, b2r20, xlim=x, main="Particles Transmitted over Time, b=1e-5", 
     xlab="Time", ylab="Particles Transmitted (Beta = 0.15)", pch=1, cex=0.5)

plot(t3r20, b3r20, xlim=x, main="Particles Transmitted over Time, b=1e-6", 
     xlab="Time", ylab="Particles Transmitted (Beta = 0.18)", pch=1, cex=0.5)

plot(t4r20, b4r20, xlim=x, main="Particles Transmitted over Time, b=1e-7", 
     xlab="Time", ylab="Particles Transmitted (Beta = 0.23)", pch=1, cex=0.5)

plot(t5r20, b5r20, xlim=x, main="Particles Transmitted over Time, b=1e-8", 
     xlab="Time", ylab="Particles Transmitted (Beta = 0.32)", pch=1, cex=0.5)

plot(t6r20, b6r20, xlim=x, main="Particles Transmitted over Time, b=1e-9", 
     xlab="Time", ylab="Particles Transmitted (Beta = 0.50)", pch=1, cex=0.5)

plot(t7r20, b7r20, xlim=x, main="Particles Transmitted over Time, b=1e-10", 
     xlab="Time", ylab="Particles Transmitted (Beta = 1.30)", pch=1, cex=0.5)



b <- c(1e-10, 1e-9, 1e-8, 1e-7, 1e-6)
TheoR <- c(2, 2, 2, 2, 2)
plot(log10(b), TheoR, ylim=c(2.0,2.2), main="Theoretical R0 vs Simulation", 
     xlab="Log10(b)", ylab="R0", pch=3, cex=0.75)

points(log10(b), R0, type="p", col="red")
lines(log10(b), R0, lty=22, col="red")
lines(log10(b), TheoR, lty=1, col="black")
legend(-10,2.20, c("Simulated", "Theoretical"), # puts text in the legend 
       
       lty=c(22, 1), cex=0.6, # gives the legend appropriate symbols (lines)
       
       col=c("red", "black")) # gives the legend lines the correct color and width