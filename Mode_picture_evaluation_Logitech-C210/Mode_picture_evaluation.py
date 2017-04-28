# Coded by AO

# +-+-+ GENERAL REMARKS +-+-+
# ---------------------------
# This program evaluates pictures taken with a Logitech C210 CCD camera.
# The edge length of one pixel is 3.3 micron.
# The program calculates the center of mass by adding all the lines in both directions
# - x and y - and fitting a Gaussian function in order to get the maximum in both
# directions and the standard deviation, which determines the with of the measured beam.
# The program outputs 5 files:
#   - Fitparameters-#-####-##-##.txt: Containing the parameters of the Gaussian functions
#   - Plot_(mm)-#-####-##-##.pdf: Plot of the picture with mm labels and widths of the fits
#   - Plot_(px)-#-####-##-##.pdf: Plot of the picture with px labels and width of the fits
#   - Xfit-#-####-##-##.pdf: Plot of the Gaussion function where all lines in y-direction
#                               were added
#   - Yfit-#-####-##-##.pdf: Plot of the Gaussian function where all lines in x-direction
#                               were added
# There are two command line options that can be deployed as standard command line program options
# with a '-' in front. The option '--source' can specify an external hard drive as a source for
# mode pictures. The option '--folder' specifies the actual folder containing the data.


import numpy as np
import matplotlib.pyplot as plt
import scipy.ndimage as sc
import scipy.signal as sig
import scipy.optimize as opt
import math
import os
import sys
import platform 
from matplotlib.ticker import FuncFormatter
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import MultipleLocator
from optparse import OptionParser

def f(x,x0,s,A,o):
	return o+(A/(2*math.pi*s**2))*np.exp(-(x-x0)**2/(2*s**2))

def cm2inch(value):
    return value/2.54

def formatfunc(x,pos):
	return "%.3f" %(x*3.3e-3)

def formatfunc640(x,pos):
    return "%.3f" %((x-300)*3.3e-3)

def formatfunc480(x,pos):
    return "%.3f" %((x-200)*3.3e-3)



parser = OptionParser()
parser.add_option("-s", "--source", action="store", dest="source")
parser.add_option("-f", "--folder", action="store", dest="folder")
(options,args) = parser.parse_args()


os.system("clear")
print "\n\n\n"
print "-------------------------------------------"
print "------------ PICTURE EVALUATION -----------"
print "-------------------------------------------"
print "\n\n\n"


foldername = options.folder+"/"

# This part needs to be adjusted for different computers or data storage systems
if (options.source == "LABDATA"):
    if (platform.system() == "Darwin"):
        path = "/Volumes/LABDATA/2016/Mode_spectra_measurements/"

    elif (platform.system() == "Linux"):
        path = ""
else:
    if (platform.system() == "Darwin"):
        path = "/Users/aotterp/PowerFolders/Measurements/Mode_picture_measurements/"

    elif (platform.system() == "Linux"):
        path = "/home/labadmin/Documents/Alex/Measurements/Mode_picture_measurements/"

files = []
files += [each for each in os.listdir(path+foldername) if each.endswith(".jpg")]


# The for loop iterates through all jpg files within the specified folder.
# It skips files that are already processed.
for ite in range(len(files)):

    filename = files[ite]
	
    if ("Fitparameters-"+filename.split(".",2)[0][5:]+".txt" in os.listdir(path+foldername)):
        print " + File "+filename+" already processed."
        continue
    else:
        im = sc.imread(path+foldername+filename,mode="L")
            
        x = np.zeros(640,dtype="float")
        y = np.zeros(480,dtype="float")

        x0 = np.arange(0,640,1)
        y0 = np.arange(0,480,1)


        # These for loops add all the picture lines in both directions
        for ite in range(len(im)):
            x += im[ite]
        for ite in range(len(im.T)):
            y += im.T[ite]


        # The fitting is done in the following. The initial parameters might have to be
        # adjusted and are arbitrarily chosen.
        xleft = 0
        xright = 600
        x0fit = np.argmax(x)
        sxfit = 100
        Axfit = 1000000000
        oxfit = np.amin(x)
        poptx, pcovx = opt.curve_fit(f,x0[xleft:xright],x[xleft:xright],[x0fit,sxfit,Axfit,oxfit])
            
        x0fit = poptx[0]
        sxfit = poptx[1] 
        Axfit = poptx[2]
        oxfit = poptx[3] 


        yleft = 0
        yright = 450
        y0fit = np.argmax(y)
        syfit = 100
        Ayfit = 1000000000
        oyfit = np.amin(y)
        popty, pcovy = opt.curve_fit(f,y0[yleft:yright],y[yleft:yright],[y0fit,syfit,Ayfit,oyfit])
            
        y0fit = popty[0]
        syfit = popty[1] 
        Ayfit = popty[2]
        oyfit = popty[3] 	
                
            
        # The following code produces the first plot of the whole picture with px labels.
        # The usetex, etc. commands might have to be commented on Windows machines.
        plt.rc("text", usetex=True)
        plt.rc("font", **{"family":"sans-serif","sans-serif":["Helvetica"],"size":11})
        plt.rcParams["text.latex.preamble"]=["\\usepackage{siunitx}","\\usepackage[helvet]{sfmath}","\\sisetup{math-rm=\mathsf,text-rm=\sffamily}"]
        plt.rcParams["legend.fontsize"]=11
        fig = plt.figure()

        fig.set_size_inches(cm2inch(13),cm2inch(8.67))
        ax = plt.gca()
        plt.imshow(im,origin="lower")
        ax.plot(x0fit,y0fit,marker="x",markersize=10,markeredgewidth=1,color="white")
        ax.plot((x0fit-sxfit,x0fit+sxfit),(y0fit,y0fit),color="white",linewidth=1.5)
        ax.plot((x0fit,x0fit),(y0fit-syfit,y0fit+syfit),color="white",linewidth=1.5)
        ax.set_xlim(0,640)
        ax.set_ylim(0,480)

        ax.set_xlabel("Pixels in x-direction")
        ax.set_ylabel("Pixels in y-direction")
        ax.annotate("", xy=(1.05,1), xycoords="axes fraction", xytext=(1.05,0), textcoords="axes fraction", arrowprops=dict(arrowstyle="->",connectionstyle="arc3"))
        ax.annotate(r"Out-coupling angle $\Phi_\text{out}$", xy=(1.1,0.5), xycoords="axes fraction", rotation="270", va="center", ha="center")
        box = ax.get_position()
        ax.set_position([box.x0-0.025,box.y0+0.05,box.width,box.height])
        plt.savefig(path+foldername+"Plot_(px)-"+filename.split(".",2)[0][5:]+".pdf", format="pdf")
        plt.clf()
        plt.close()


        # The following code produces the second plot of the whole picture with mm labels.
        # The usetex, etc. commands might have to be commented on Windows machines.
        plt.rc("text", usetex=True)
        plt.rc("font", **{"family":"sans-serif","sans-serif":["Helvetica"],"size":11})
        plt.rcParams["text.latex.preamble"]=["\\usepackage{siunitx}","\\usepackage[helvet]{sfmath}","\\sisetup{math-rm=\mathsf,text-rm=\sffamily}"]
        plt.rcParams["legend.fontsize"]=11
        fig = plt.figure()

        fig.set_size_inches(cm2inch(13),cm2inch(8.67))
        ax = plt.gca()
        plt.imshow(im,origin="lower")
        ax.plot(x0fit,y0fit,marker="x",markersize=10,markeredgewidth=1,color="white")
        ax.plot((x0fit-sxfit,x0fit+sxfit),(y0fit,y0fit),color="white",linewidth=1.5)
        ax.plot((x0fit,x0fit),(y0fit-syfit,y0fit+syfit),color="white",linewidth=1.5)
        ax.set_xlim(0,640)
        ax.set_ylim(0,480)
        ax.xaxis.set_major_formatter(FuncFormatter(formatfunc))
        ax.yaxis.set_major_formatter(FuncFormatter(formatfunc))

        ax.set_xlabel(r"Length in x-direction $[\SI{}{\milli\metre}]$")
        ax.set_ylabel(r"Length in y-direction $[\SI{}{\milli\metre}]$")
        ax.annotate("", xy=(1.05,1), xycoords="axes fraction", xytext=(1.05,0), textcoords="axes fraction", arrowprops=dict(arrowstyle="->",connectionstyle="arc3"))
        ax.annotate(r"Out-coupling angle $\Phi_\text{out}$", xy=(1.1,0.5), xycoords="axes fraction", rotation="270", va="center", ha="center")
        box = ax.get_position()
        ax.set_position([box.x0-0.025,box.y0+0.05,box.width,box.height])
        plt.savefig(path+foldername+"Plot_(mm)-"+filename.split(".",2)[0][5:]+".pdf", format="pdf")
        plt.clf()
        plt.close()


        # The following code produces the Xfit plot containing the Gaussian function
        # fitted to the added lines in y-direction.
        # The usetex, etc. commands might have to be commented on Windows machines.
        plt.rc("text", usetex=True)
        plt.rc("font", **{"family":"sans-serif","sans-serif":["Helvetica"],"size":11})
        plt.rcParams["text.latex.preamble"]=["\\usepackage{siunitx}","\\usepackage[helvet]{sfmath}","\\sisetup{math-rm=\mathsf,text-rm=\sffamily}"]
        plt.rcParams["legend.fontsize"]=11
        fig = plt.figure(tight_layout=True)
        fig.set_size_inches(cm2inch(13),cm2inch(8.67))
        ax = fig.add_subplot(111)
        ax.plot(x0,x)
        ax.plot(x0,f(x0,x0fit,sxfit,Axfit,oxfit),color="red")
        ax.xaxis.set_major_formatter(FuncFormatter(formatfunc640))
        ax.set_yticks(np.delete(ax.get_yticks(),0))
        ax.set_xlabel(r"Length in x-direction $[\SI{}{\milli\metre}]$")
        ax.set_ylabel("Added lines in y-direction")
        ax.grid(True)
        ax.set_xlim(0,640)
        ax.ticklabel_format(style="sci", axis="y", scilimits=(0,0))
        ax.annotate(r"$x_0 =\SI{%.3f}{\milli\metre}$"%(poptx[0]*3.3e-3)+"\n"+r"$\sigma = \SI{%.3f}{\milli\metre}$"%(poptx[1]*3.3e-3),xy=(1,1),xycoords="axes fraction", xytext=(-10,-10), textcoords="offset points", ha="right", va="top", bbox=dict(boxstyle="square,pad=0.3",ec="black",fc="white"))
        plt.savefig(path+foldername+"Xfit-"+filename.split(".",2)[0][5:]+".pdf", format="pdf")
        plt.clf()
        plt.close()


        # The following code produces the Yfit plot containing the Gaussian function
        # fitted to the added lines in x-direction.
        # The usetex, etc. commands might have to be commented on Windows machines.
        plt.rc("text", usetex=True)
        plt.rc("font", **{"family":"sans-serif","sans-serif":["Helvetica"],"size":11})
        plt.rcParams["text.latex.preamble"]=["\\usepackage{siunitx}","\\usepackage[helvet]{sfmath}","\\sisetup{math-rm=\mathsf,text-rm=\sffamily}"]
        plt.rcParams["legend.fontsize"]=11
        fig = plt.figure(tight_layout=True)
        fig.set_size_inches(cm2inch(13),cm2inch(8.67))
        ax = fig.add_subplot(111)
        ax.plot(y0,y)
        ax.plot(y0,f(y0,y0fit,syfit,Ayfit,oyfit),color="red")
        ax.xaxis.set_major_formatter(FuncFormatter(formatfunc480))
        ax.set_xlabel(r"Length in y-direction $[\SI{}{\milli\metre}]$")
        ax.set_ylabel("Added lines in x-direction")
        ax.grid(True)
        ax.set_xlim(0,480)
        ax.annotate("", xy=(0.2,1.05), xycoords="axes fraction", xytext=(0.8,1.05), textcoords="axes fraction", arrowprops=dict(arrowstyle="<-",connectionstyle="arc3"))
        ax.annotate(r"$\Phi_\text{out}$", xy=(0.84,1.05), xycoords="axes fraction", va="center", ha="center")
        ax.ticklabel_format(style="sci", axis="y", scilimits=(0,0))
        ax.annotate(r"$y_0 =\SI{%.3f}{\milli\metre}$"%(popty[0]*3.3e-3)+"\n"+r"$\sigma = \SI{%.3f}{\milli\metre}$"%(popty[1]*3.3e-3),xy=(1,1),xycoords="axes fraction", xytext=(-10,-10), textcoords="offset points", ha="right", va="top", bbox=dict(boxstyle="square,pad=0.3",ec="black",fc="white"))
        plt.savefig(path+foldername+"Yfit-"+filename.split(".",2)[0][5:]+".pdf", format="pdf")
        plt.close()


        # The following code stores all fitting parameters in a separate file. 
        outputfile = open(path+foldername+"Fitparameters-"+filename.split(".",2)[0][5:]+".txt", "w")
        outputfile.write("# Xfit Parameters x0fit, sxfit, Axfit and oxfit:\n")
        outputfile.write(str(x0fit)+"\n")
        outputfile.write(str(sxfit)+"\n")
        outputfile.write(str(Axfit)+"\n")
        outputfile.write(str(oxfit)+"\n")
        outputfile.write("# Yfit Parameters y0fit, syfit, Ayfit and oyfit:\n")
        outputfile.write(str(y0fit)+"\n")
        outputfile.write(str(syfit)+"\n")
        outputfile.write(str(Ayfit)+"\n")
        outputfile.write(str(oyfit)+"\n")
        outputfile.close()

        print " + File "+filename+" processed."

    
    








