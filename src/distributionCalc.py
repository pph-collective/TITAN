from scipy.stats import norm,rayleigh
from numpy import linspace
from pylab import plot,show,hist,figure,title

samp = rayleigh.rvs(loc=5,scale=2,size=10) # samples generation
print samp
samp = [1,5,5,5,5,11,11,11,11,11,23,23,24]

param = rayleigh.fit(samp) # distribution fitting

x = linspace(5,13,100)
# fitted distribution
pdf_fitted = rayleigh.pdf(x,loc=param[0],scale=param[1])
# original distribution
pdf = rayleigh.pdf(x,loc=5,scale=2)

title('Rayleigh distribution')
plot(x,pdf_fitted,'r-',x,pdf,'b-')
hist(samp,normed=1,alpha=.3)
show()