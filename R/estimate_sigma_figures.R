# plot dispersal curves

s = 0.308 # estimated sigma 
sll = 0.187 # lower 95% CI
sul = 0.441 # upper 95% CI
x = seq(0,2,0.01)


px = function(x,s, verbose = FALSE){(1/(s*sqrt(2)))*exp(-sqrt(2)*x/s)} # # laplacian kernal
s_lapKernal = px(x,s); s_lapKernal_normalized = s_lapKernal/max(s_lapKernal)
sll_lapKernal = px(x,sll); sll_lapKernal_normalized = sll_lapKernal/max(sll_lapKernal)
sul_lapKernal = px(x,sul); sul_lapKernal_normalized = sul_lapKernal/max(sul_lapKernal)
plot(x,s_lapKernal_normalized,type="l",ylab="Dispersal strength (normalized)",xlab="Distance (km)")
points(x,sll_lapKernal_normalized,type="l",col="red")
points(x,sul_lapKernal_normalized,type="l",col="red")
print(paste("Mean dispersal with Laplacian = ",s/sqrt(2)))
segments(x0 = s/sqrt(2),y0=-1,x1=s/sqrt(2),y1=0.5,lty="dashed")



