xlabR0tot="Simulated R0"
xlabR0V1= "Simulated R0 for Type 1"
xlabR0V2="Simulated R0"
ylabMW="Mean Bottleneck Width"
plot(MWvR0, xlab=xlabR0V2, ylab=ylabMW, type = "o", col="blue")
lines(MWvR01, type="o", col="red")
title(main="Bottleneck Width as a Function of Simulated R0", col.main="Black", font.main=1, cex.main=0.9, cex.sub=0.65, sub="Theoretical single strain type 1 R0 fixed at 0.60, type 2 R0 fixed at 1.40")
legend(0.15,126, c("Type 1", "Type 2"), # puts text in the legend 
       
       lty=c(1, 1), cex=0.6, # gives the legend appropriate symbols (lines)
       
       col=c("blue", "red")) # gives the legend lines the correct color and width
