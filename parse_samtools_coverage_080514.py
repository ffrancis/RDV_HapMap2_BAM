# Parse samtool mpileup output to calculate coverage per specified bins
# Change x in (kmer_size = x) to set the desired bin width.
# Make sure there are no duplicates in the mpileup output (with respect to coordinates)

#import re
import matplotlib.pyplot as plt
import pandas as pd
import math

#filename1 = "testcov_CML322_chr3_1-10.txt"
#filename2 = "testcov_B73_chr3_1-10.txt"

filename1 = "cov_CML322_chr3_219833000-219905001.txt"                                      #samtool mpileup output
filename2 = "cov_B73_chr3_219833000-219905001.txt"  


##########
#filename1
##########
filename_readname=filename1.replace('-','_').replace('.','_')
reg_start = filename_readname.split("_")[3]                                  #get region start and stop coordinates from filename
reg_stop = filename_readname.split("_")[4]
#print str(reg_start) + "\t" + str(reg_stop)

f1 = open(filename1, "r")
lines1 = f1.readlines()

coord = []
no_reads= []
for line in lines1:                                                  # get coordinates and respective number of reads from cov file from samtool mpileup                                                     
    line2=line.split('\n')
    #print line2 [0]
    line3 = line2 [0]
    cols = line3.split("\t")
    if ((int(cols[1]))>= int(reg_start)) and ((int(cols[1]))<= int(reg_stop)):
        no_reads.append(int(cols[3]))
        coord.append(int(cols[1]))
#print (coord)
#print (no_reads)

df_cov = pd.DataFrame(index=coord, data=no_reads )                    #create a pandas data frame(df_cov) for coverage
#print df_cov
ref_coord = range(int(reg_start), int(reg_stop)+1)
#print (ref_coord)
column = [int(0) for number in xrange(int(reg_stop)-int(reg_start)+1)]
#print column
df_ref_coord = pd.DataFrame(index = ref_coord, data =column)                       #create a pandas data frame(df_ref_coord) (only index) for b73 ref coordinates based on start stop position from file name
df_merged = (df_ref_coord +df_cov)                                                 #add two dataframes(only possible if both df's have same index and colum names)
df_merged = df_merged.fillna(0)                                                    # change NaN to 0
df_merged.columns = ['no_reads']

#print list(df_merged.no_reads)

x= (list(df_merged.no_reads))
#print x
f1.close()

##########
#filename2
##########
filename_readname=filename2.replace('-','_').replace('.','_')
reg_start = filename_readname.split("_")[3]                                  #get region start and stop coordinates from filename
reg_stop = filename_readname.split("_")[4]
#print str(reg_start) + "\t" + str(reg_stop)

f1 = open(filename2, "r")
lines1 = f1.readlines()

coord = []
no_reads= []
for line in lines1:                                                  # get coordinates and respective number of reads from cov file from samtool mpileup                                                     
    line2=line.split('\n')
    #print line2 [0]
    line3 = line2 [0]
    cols = line3.split("\t")
    if ((int(cols[1]))>= int(reg_start)) and ((int(cols[1]))<= int(reg_stop)):
        no_reads.append(int(cols[3]))
        coord.append(int(cols[1]))
#print (coord)
#print (no_reads)

df_cov = pd.DataFrame(index=coord, data=no_reads )                    #create a pandas data frame(df_cov) for coverage
#print df_cov
ref_coord = range(int(reg_start), int(reg_stop)+1)
#print (ref_coord)
column = [int(0) for number in xrange(int(reg_stop)-int(reg_start)+1)]
#print column
df_ref_coord = pd.DataFrame(index = ref_coord, data =column)                       #create a pandas data frame(df_ref_coord) (only index) for b73 ref coordinates based on start stop position from file name
df_merged = (df_ref_coord +df_cov)                                                 #add two dataframes(only possible if both df's have same index and colum names)
df_merged = df_merged.fillna(0)                                                    # change NaN to 0
df_merged.columns = ['no_reads']

y= (list(df_merged.no_reads))
#print len(y)
#print len(x)
f1.close()

##########
#PLOT -simple scatter
##########
#log2_ratio = [(math.log((a+0.5),2))/(math.log((b+0.5),2)) for a, b in zip(x, y)]     #log of zero gives us an infinite number, negative infinity. So one trick is to add a small value such as a half. And then we get negative 1 instead of negative infinity. So it helps us to plot these zero values.(only for plotting)
log2_ratio = [(math.log(((a+.5)/(b+.5)),2)) for a, b in zip(x, y)]

#log2_ratio = map(int, log2_ratio)
#ref_coord = map(int, ref_coord)

print "log2ratio"+"\t"+ str(len(log2_ratio))
print "ref_coord"+"\t"+str(len(ref_coord))

#plt.scatter(ref_coord, log2_ratio)
##plt.plot([int(reg_start), int(reg_stop)], [0, 0], 'k-', lw=2)
#plt.show()


##########
##Hexbin PLot
##########
from astroML.plotting import setup_text_plots
setup_text_plots(fontsize=20, usetex=True)
plt.hexbin(ref_coord, log2_ratio, cmap=plt.cm.binary, gridsize=48, norm=None)       #72001/48 ~= 1500 bp less than the average size of a gene in maize; change this for each region
cb = plt.colorbar()
plt.title("Read depth plot CML322 vs B73 chr3:219833000-219905001")
plt.xlabel(r'${\rm Chromosome - position}$')
plt.ylabel(r'${\rm Log_2 - Ratio}$')
plt.show()




###########
##PLOT- contour scatter
##http://www.astroml.org/book_figures/chapter1/fig_S82_scatter_contour.html
###########
#from astroML.plotting import scatter_contour
#from astroML.datasets import fetch_sdss_S82standards
#from astroML.plotting import setup_text_plots
#setup_text_plots(fontsize=20, usetex=True)
#
#from numpy  import *                                                                   #Change list to array; for contour plotting
#a = array(ref_coord)
#b = array(log2_ratio)
#
##print type(ref_coord )
##print type(log2_ratio)
#print type(a)
#print type(b)
#
### plot the results
#fig, ax = plt.subplots(figsize=(15, 13.75))
#scatter_contour(a, b, threshold=100, log_counts=True, ax=ax,
#                histogram2d_args=dict(bins=48),                                        #72001/48 ~= 1500 bp less than the average size of a gene in maize; change this for each region
#                plot_args=dict(marker=',', linestyle='none', color='black'),
#                contour_args=dict(cmap=plt.cm.bone))
#
#
#ax.set_xlabel(r'${\rm Chromosome - position}$')
#ax.set_ylabel(r'${\rm Log_2 - Ratio}$')
#
#ax.set_title("Read depth plot CML322 vs B73 chr3:219833000-219905001")                                      #Change for each region                                                                         
#ax.set_xlim(219833000, 219905001)                                                                           #Change for each region  
##ax.set_ylim(-0.6, 2.5)
#plt.show()

