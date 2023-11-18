import numpy as np
from scipy.stats import expon
from scipy.optimize import curve_fit

# load data
import csv
with open("peaks_cover_sample_count_ns_2018.csv",'r') as inf:
	lines = inf.readlines()
	y = []
	for line in lines:
		y.append(int(line.split()[1]))
x = np.array([i for i in range(1,len(y)+1)])

# define the exponential function
def exp_func(x, a, b, c):
    return a * np.exp(-b * (x+c))

# fit the data to the exponential function
params, cov = curve_fit(exp_func, x[1:21], y[1:21], maxfev=200000)
print(params)

# get the fitted values
fitted_values = exp_func(x, *params)

# plot the data and the fitted curve
import matplotlib
#set font
matplotlib.rcParams["font.sans-serif"] = ["Arial"]

import matplotlib.pyplot as plt
# set the plot size
plt.figure(figsize=(4, 4))
# set plot range
plt.xlim(0,len(y))
# set log scale
plt.yscale('log')
# set xy label
plt.xlabel('Occupancy score',fontsize=12)
plt.ylabel('Number of origin sites',fontsize=12)
# plot the data points
plt.scatter(x,y, label='Observed', color='navy',s = 4, alpha = 0.8)
# plot the fitted curve
plt.plot(fitted_values, label='Background',color='limegreen')
# set title
plt.title("Exponential model fitting\nfor background noise", fontsize=14)
# plot legend
plt.legend()
# plot cutoff line
plt.axvline(x=20,ymax=1000,ymin=0,color='grey',linestyle='--')
# Add a label next to the abline
plt.text(20, 100, 'Occupancy score=20\nFDR=0.1', color='grey')
# save plot
plt.savefig('origins_occupancy_score_exp_model_fit_ns_2-21.pdf',bbox_inches = 'tight',pad_inches = 0.1,dpi=600,transparent=True)
plt.show()
plt.close()

## print out the oddsRatio
sum_ori,sum_fit,flag = 0,0,0
yfit = fitted_values
for i in np.arange(len(x)-1,1,-1):
    sum_ori += y[i]
    sum_fit += yfit[i]
    #flag = max(flag,sum_ori/sum_fit)
    print(i,y[i],yfit[i],(sum_ori)/sum_fit)
