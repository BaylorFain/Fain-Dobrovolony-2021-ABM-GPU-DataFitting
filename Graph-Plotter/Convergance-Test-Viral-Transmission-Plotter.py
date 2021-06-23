#This program plots the viral titer curves for each individual run, the median curves of all the runs, the histograms for each measurable aspect, the measurable aspect, and every run and the median curves on one plot. It also calculates if the distribution of the parameters is normal.

import numpy
import datetime
import matplotlib
import matplotlib.pyplot as plt
import os
from scipy.stats import norm
from scipy.stats import shapiro
from scipy.stats import anderson
from scipy.stats import normaltest
from scipy.stats import ttest_1samp
import scipy.integrate as integrate
import matplotlib.ticker as ticker
from matplotlib.ticker import FormatStrFormatter

os.system('cls' if os.name=='nt' else 'clear')

def growthfit(TimeArray, VirusArray, IndexOfPeakTime):
    import scipy.optimize as optim
    import numpy as np
    if IndexOfPeakTime <= 2:
        b = 0.0
        perr = [0.0, 0.0, 0.0] 
        r2 = 0.0
    else:
        def logistic(t, a, b, c):
            return a/(1+np.exp(-b*(t-c)))
        beginningtime = int(numpy.nanargmin(VirusArray[0:IndexOfPeakTime])+6)
        x = TimeArray[beginningtime:IndexOfPeakTime]
        y = np.divide(VirusArray[beginningtime:IndexOfPeakTime],VirusArray[IndexOfPeakTime])
        bounds = (0, [2,1,IndexOfPeakTime])
        (a,b,c),cov = optim.curve_fit(logistic, x, y, bounds=bounds)
        perr = np.sqrt(np.diag(cov))
        y_fit = numpy.empty([len(x)],dtype=float )
        for m in range(len(x)):
            y_fit[m] = logistic(x[m], a, b, c)
        # residual sum of squares
        ss_res = np.sum((y - y_fit) ** 2)
        # total sum of squares
        ss_tot = np.sum((y - np.mean(y)) ** 2)
        # r-squared
        r2 = 1 - (ss_res / ss_tot)
    return (b, perr[1], r2)

def decayfit(TimeArray,EclipseArray,InfectedArray,DeadArray,VirusArray):
    import scipy.optimize as optim
    import numpy as np
    notif = 0
    StartTime = 0
    EndTime = TimeArray[-1]
    for ii in range(len(TimeArray)):
        if (EclipseArray[ii]==0)&(InfectedArray[ii]==0)&(DeadArray[ii]!=0)&(notif==0):
            StartTime = TimeArray[ii]
            notif = 1
        if (VirusArray[ii] <= 10**2)&(TimeArray[ii] > 5):
            EndTime = TimeArray[ii]
            break

    x = TimeArray[StartTime:EndTime]
    y = np.log(VirusArray[StartTime:EndTime])    

    def linear(x, m, b):
        return m*x+b

    try:
        (m,b),cov = optim.curve_fit(linear, x, y)
        perr = np.sqrt(np.diag(cov))

        y_fit = numpy.empty([len(x)],dtype=float )
        for iii in range(len(x)):
            y_fit[iii] = linear(iii,m,b)
        # residual sum of squares
        ss_res = np.sum((y - y_fit) ** 2)
        # total sum of squares
        ss_tot = np.sum((y - np.mean(y)) ** 2)
        # r-squared
        r2 = 1 - (ss_res / ss_tot)
        return(np.abs(m), perr[1], r2)
    except ValueError:
        return([np.nan])

MOI = [numpy.power(10.0,-2)]
logMOI = numpy.round(numpy.log10(MOI),1)

TimeSteps = numpy.array(["0.000625", "0.001280", "0.002500", "0.005000", "0.010000", "0.020000", "0.040000"])

NumberOfRuns = 10

#Cellfree-InitialCell
#OuterPath = r"/home/baylor/Documents/Research/Papers/Current_Papers/CompSci-Paper/Graph-Plotter/Convergance/Cellfree-InitialCell/Cellfree-InitialCell-Runs"

#Cellfree-InitialVirus
#OuterPath = r"/home/baylor/Documents/Research/Papers/Current_Papers/CompSci-Paper/Graph-Plotter/Convergance/Cellfree-InitialVirus/Cellfree-InitialVirus-Runs"

#Neither-InitialVirus
OuterPath = r"/home/baylor/Documents/Research/Papers/Current_Papers/CompSci-Paper/Graph-Plotter/Convergance/Neither-InitialVirus/Neither-InitialVirus-Runs"

BigPeakVirus = numpy.empty([len(TimeSteps), NumberOfRuns],dtype=float )
BigTimeOfPeak = numpy.empty([len(TimeSteps), NumberOfRuns],dtype=float )
BigUpSlope = numpy.empty([len(TimeSteps), NumberOfRuns],dtype=float )
BigDownSlope = numpy.empty([len(TimeSteps), NumberOfRuns],dtype=float )
BigEndTime = numpy.empty([len(TimeSteps), NumberOfRuns],dtype=float )
BigAUC = numpy.empty([len(TimeSteps), NumberOfRuns],dtype=float )
        
BigTimeArray = numpy.empty([NumberOfRuns],dtype=tuple )
BigHealthyArray = numpy.empty([NumberOfRuns],dtype=tuple )
BigEclipseArray = numpy.empty([NumberOfRuns],dtype=tuple )
BigInfectedArray = numpy.empty([NumberOfRuns],dtype=tuple )
BigDeadArray = numpy.empty([NumberOfRuns],dtype=tuple )
BigVirusArray = numpy.empty([NumberOfRuns],dtype=tuple )

ReallyBigVirusArray = numpy.empty([len(TimeSteps),NumberOfRuns],dtype=tuple )
BigMedianVirusArray = numpy.empty([len(TimeSteps)],dtype=tuple )

LongestTime = 0.0

for j in range(len(TimeSteps)):
    OUTERPATH = os.path.join(OuterPath).replace('\\','//')
    directory = os.path.join("".join(OUTERPATH+ "_" +str(TimeSteps[j])+ "_" +"Analysis"))
    if not os.path.exists(directory):
        os.makedirs(directory)
    OUTER_Path_to_Folder = os.path.abspath(directory)
    
    directory = os.path.join("".join(OUTERPATH+ "_" +"Graphs"))
    if not os.path.exists(directory):
        os.makedirs(directory)
    OUTEROUTER_Path_to_Folder = os.path.abspath(directory)
    
    for i in range(NumberOfRuns):    
        os.system('cls' if os.name=='nt' else 'clear')
        print(str(TimeSteps[j])+"\n")
        print(i)
        
#        CELLFREE
#        InnerPath = r"/607_"+str(i)+"-CELLFREE_"+str(logMOI[0])+"-MOI_"+str(TimeSteps[j])+"-timestep"

#        NEITHER
        InnerPath = r"/607_"+str(i)+"-Neither_"+str(logMOI[0])+"-MOI_"+str(TimeSteps[j])+"-timestep"
        
        INNERPATH = os.path.join(OuterPath+InnerPath).replace('\\','//')
        directory = os.path.join("".join(INNERPATH + "/Analysis"))
        if not os.path.exists(directory):
            os.makedirs(directory)
        Inner_Path_to_Folder = os.path.abspath(directory)

        TimeArray = []
        HealthyArray = []
        EclipseArray = []
        InfectedArray = []
        DeadArray = []
        VirusArray = []

        with open(os.path.join(INNERPATH,"PerTimeStep.txt")) as CellInFile:
            for cellline in CellInFile:
                CellData = cellline
                CellData = CellData.split(",")
                del CellData[-1]
                
                if(LongestTime < int(CellData[0])):
                    LongestTime = int(CellData[0])
                    
                TimeArray.append(int(CellData[0]))
                HealthyArray.append(int(CellData[1]))
                EclipseArray.append(int(CellData[2]))
                InfectedArray.append(int(CellData[3]))
                DeadArray.append(int(CellData[4]))
                VirusArray.append(float(CellData[5]))

        BigTimeArray[i] = TimeArray
        BigHealthyArray[i] = HealthyArray
        BigEclipseArray[i] = EclipseArray
        BigInfectedArray[i] = InfectedArray
        BigDeadArray[i] = DeadArray
        BigVirusArray[i] = VirusArray

        ReallyBigVirusArray[j][i] = VirusArray

##################################################################################
##Plots each run in the analysis folder in that runs folder. Takes up the most time
##################################################################################
#        plt.rcParams['figure.figsize'] = [10, 10]
#        fig, ax = plt.subplots()
#        plt.plot(TimeArray, HealthyArray, linewidth = 3)
#        plt.xticks(numpy.arange(0, TimeArray[-1], (TimeArray[-1]/10)))
#        plt.tight_layout()
#        plt.close()
#        fig.savefig(os.path.join(Inner_Path_to_Folder,"HealthyVsTime"+".pdf"), dpi=fig.dpi, bbox_inches='tight')
#                    
#        plt.rcParams['figure.figsize'] = [10, 10]
#        fig, ax = plt.subplots()
#        plt.plot(TimeArray, EclipseArray, linewidth = 3)
#        plt.xticks(numpy.arange(0, TimeArray[-1], (TimeArray[-1]/10)))
#        plt.tight_layout()
#        plt.close()
#        fig.savefig(os.path.join(Inner_Path_to_Folder,"EclipseVsTime"+".pdf"), dpi=fig.dpi, bbox_inches='tight')
#                    
#        plt.rcParams['figure.figsize'] = [10, 10]
#        fig, ax = plt.subplots()
#        plt.plot(TimeArray, InfectedArray, linewidth = 3)
#        plt.xticks(numpy.arange(0, TimeArray[-1], (TimeArray[-1]/10)))
#        plt.tight_layout()
#        plt.close()
#        fig.savefig(os.path.join(Inner_Path_to_Folder,"InfectedVsTime"+".pdf"), dpi=fig.dpi, bbox_inches='tight')

#        plt.rcParams['figure.figsize'] = [10, 10]
#        fig, ax = plt.subplots()
#        plt.plot(TimeArray, DeadArray, linewidth = 3)
#        plt.xticks(numpy.arange(0, TimeArray[-1], (TimeArray[-1]/10)))
#        plt.close()
#        fig.savefig(os.path.join(Inner_Path_to_Folder,"DeadVsTime"+".pdf"), dpi=fig.dpi, bbox_inches='tight')

#        plt.rcParams['figure.figsize'] = [10, 10]
#        fig, ax = plt.subplots()
#        plt.yscale("log")
#        plt.plot(TimeArray, VirusArray, linewidth = 3)
#        plt.xticks(numpy.arange(0, TimeArray[-1], (TimeArray[-1]/10)))
#        plt.tight_layout()
#        plt.title("Virus vs Time")
#        plt.ylabel("Amount of Virus")
#        plt.xlabel("Time-(hr)")
#        plt.close()
#        fig.savefig(os.path.join(Inner_Path_to_Folder,"VirusVsTime"+".pdf"), dpi=fig.dpi, bbox_inches='tight')

##################################################################################        
        BeginingTime = [[-1]]#numpy.argwhere(numpy.array(BigVirusArray[i]) == 0.0)#
        IndexOfPeakTime = numpy.nanargmax(BigVirusArray[i])

        if j == 6:
            BigPeakVirus[j][i] = numpy.nan
            BigTimeOfPeak[j][i] = numpy.nan
            BigUpSlope[j][i] = numpy.nan
            BigDownSlope[j][i] = numpy.nan
            BigEndTime[j][i] = numpy.nan
            BigAUC[j][i] = numpy.nan
        else:
            BigPeakVirus[j][i] = max(BigVirusArray[i])
            BigTimeOfPeak[j][i] = BigTimeArray[i][IndexOfPeakTime]
            BigUpSlope[j][i] = growthfit(BigTimeArray[i], BigVirusArray[i], IndexOfPeakTime)[0]
            BigDownSlope[j][i] = decayfit(BigTimeArray[i], BigEclipseArray[i], BigInfectedArray[i], BigDeadArray[i], BigVirusArray[i])[0]
            BigEndTime[j][i] = BigTimeArray[i][-1]
            BigAUC[j][i] = numpy.trapz(numpy.log(BigVirusArray[i][(BeginingTime[-1][0]+1):]))
            
    MedianTimeArray = []
    MedianHealthyArray = []
    MedianEclipseArray = []
    MedianInfectedArray = []
    MedianDeadArray = []
    MedianVirusArray = []

    InbetweenTimeArray = []
    InbetweenHealthyArray = []
    InbetweenEclipseArray = []
    InbetweenInfectedArray = []
    InbetweenDeadArray = []
    InbetweenVirusArray = []

    for i in range(LongestTime):
        for k in range(NumberOfRuns):       
            try:
                InbetweenTimeArray.append(BigTimeArray[k][i])
            except IndexError:
                InbetweenTimeArray.append(0.0)
                
            try:
                InbetweenHealthyArray.append(BigHealthyArray[k][i])
            except IndexError:
                InbetweenHealthyArray.append(0.0)
                
            try:
                InbetweenEclipseArray.append(BigEclipseArray[k][i])
            except IndexError:
                InbetweenEclipseArray.append(0.0)
                
            try:
                InbetweenInfectedArray.append(BigInfectedArray[k][i])
            except IndexError:
                InbetweenInfectedArray.append(0.0)
                
            try:
                InbetweenDeadArray.append(BigDeadArray[k][i])
            except IndexError:
                InbetweenDeadArray.append(0.0)
                
            try:
                InbetweenVirusArray.append(BigVirusArray[k][i])
            except IndexError:
                InbetweenVirusArray.append(0.0)
                
        MedianTimeArray.append(numpy.median(InbetweenTimeArray))
        MedianHealthyArray.append(numpy.median(InbetweenHealthyArray))
        MedianEclipseArray.append(numpy.median(InbetweenEclipseArray))
        MedianInfectedArray.append(numpy.median(InbetweenInfectedArray))
        MedianDeadArray.append(numpy.median(InbetweenDeadArray))
        MedianVirusArray.append(numpy.median(InbetweenVirusArray))
        
        InbetweenTimeArray.clear()
        InbetweenHealthyArray.clear()
        InbetweenEclipseArray.clear()
        InbetweenInfectedArray.clear()
        InbetweenDeadArray.clear()
        InbetweenVirusArray.clear()

    BigMedianVirusArray[j] = MedianVirusArray  

    with open(os.path.join(OUTER_Path_to_Folder,str(TimeSteps[j])+"MedianPerTimeStep.txt"), 'w') as outfile:
        for i in range(len(MedianTimeArray)):
            print(str(int(MedianTimeArray[i]))+", "+str(int(MedianHealthyArray[i]))+", "+str(int(MedianEclipseArray[i]))+", "+str(int(MedianInfectedArray[i]))+", "+str(int(MedianDeadArray[i]))+", "+str(MedianVirusArray[i])+",",file=outfile)

##################################################################################
##Plots the median curve with the individual curves
##################################################################################    
#    plt.rcParams['figure.figsize'] = [10, 10]
#    fig, ax = plt.subplots()
#    for i in BigHealthyArray:
#        plt.plot(i, color = "gray", alpha = 0.3)
#    plt.plot(MedianHealthyArray, color = "red")
#    plt.xticks(numpy.arange(0, LongestTime, LongestTime/10))
#    plt.tight_layout()
#    plt.close()
#    fig.savefig(os.path.join(OUTER_Path_to_Folder,"TimeStep-" + str(TimeSteps[j]) + "_Median_HealthyVsTime"+".pdf"), dpi=fig.dpi, bbox_inches='tight')
#                
#    plt.rcParams['figure.figsize'] = [10, 10]
#    fig, ax = plt.subplots()
#    for i in BigEclipseArray:
#        plt.plot(i, color = "gray", alpha = 0.3)
#    plt.plot(MedianEclipseArray, color = "red")
#    plt.xticks(numpy.arange(0, LongestTime, LongestTime/10))
#    plt.tight_layout()
#    plt.close()
#    fig.savefig(os.path.join(OUTER_Path_to_Folder,"TimeStep-" + str(TimeSteps[j]) + "Median_EclipseVsTime"+".pdf"), dpi=fig.dpi, bbox_inches='tight')
#                
#    plt.rcParams['figure.figsize'] = [10, 10]
#    fig, ax = plt.subplots()
#    for i in BigInfectedArray:
#        plt.plot(i, color = "gray", alpha = 0.3)
#    plt.plot(MedianInfectedArray, color = "red")
#    plt.xticks(numpy.arange(0, LongestTime, LongestTime/10))
#    plt.tight_layout()
#    plt.close()
#    fig.savefig(os.path.join(OUTER_Path_to_Folder,"TimeStep-" + str(TimeSteps[j]) + "Median_InfectedVsTime"+".pdf"), dpi=fig.dpi, bbox_inches='tight')

#    plt.rcParams['figure.figsize'] = [10, 10]
#    fig, ax = plt.subplots()
#    for i in BigDeadArray:
#        plt.plot(i, color = "gray", alpha = 0.3)
#    plt.plot(MedianDeadArray, color = "red")
#    plt.xticks(numpy.arange(0, LongestTime, LongestTime/10))
#    plt.tight_layout()
#    plt.close()
#    fig.savefig(os.path.join(OUTER_Path_to_Folder,"TimeStep-" + str(TimeSteps[j]) + "Median_DeadVsTime"+".pdf"), dpi=fig.dpi, bbox_inches='tight')

#    plt.rcParams['figure.figsize'] = [10, 10]
#    fig, ax = plt.subplots()
#    plt.yscale("log")
#    for i in BigVirusArray:
#        plt.plot(i, color = "gray", alpha = 0.3)
#    plt.plot(MedianVirusArray, color = "red")
#    plt.xticks(numpy.arange(0, LongestTime, LongestTime/10))
#    plt.tight_layout()
#    plt.title("Median Virus vs Time")
#    plt.ylabel("Median_ Virus Titer")
#    plt.xlabel("Time-(hr)")
#    plt.close()
#    fig.savefig(os.path.join(OUTER_Path_to_Folder,"TimeStep-" + str(TimeSteps[j]) + "Median_VirusVsTime"+".pdf"), dpi=fig.dpi, bbox_inches='tight')

##################################################################################
##Plots a histogram for each measurable aspect of a viral titer curve
##################################################################################
##BigPeakVirus
#    histbins = numpy.linspace(min(BigPeakVirus[j]),max(BigPeakVirus[j]),20)
#    hist, bins = numpy.histogram(BigPeakVirus[j], bins=histbins, density=True)
#    width = 0.7 * (bins[1] - bins[0])
#    center = (bins[:-1] + bins[1:]) / 2
#    
#    mean, STD = norm.fit(BigPeakVirus[j])
#    x = numpy.linspace(min(BigPeakVirus[j]), max(BigPeakVirus[j]), 1000)
#    p = norm.pdf(x, mean, STD)
#    
#    plt.rcParams['figure.figsize'] = [10, 10]
#    fig, ax = plt.subplots()
#    ax.bar(center, hist, align='center', width=width)
#    plt.plot(x, p)
#    plt.tight_layout()
#    plt.close()
#    fig.savefig(os.path.join(OUTER_Path_to_Folder,str(TimeSteps[j])+"_PeakVirusHist"+".pdf"), dpi=fig.dpi, bbox_inches='tight')

##BigTimeOfPeak
#    histbins = numpy.linspace(min(BigTimeOfPeak[j]),max(BigTimeOfPeak[j]),20)
#    hist, bins = numpy.histogram(BigTimeOfPeak[j], bins=histbins, density=True)
#    width = 0.7 * (bins[1] - bins[0])
#    center = (bins[:-1] + bins[1:]) / 2
#    
#    mean, STD = norm.fit(BigTimeOfPeak[j])
#    x = numpy.linspace(min(BigTimeOfPeak[j]), max(BigTimeOfPeak[j]), 1000)
#    p = norm.pdf(x, mean, STD)
#    
#    plt.rcParams['figure.figsize'] = [10, 10]
#    fig, ax = plt.subplots()
#    ax.bar(center, hist, align='center', width=width)
#    plt.plot(x, p)
#    plt.ylabel("Frequency")
#    plt.xlabel("Time of peak virus (hours)")
#    plt.tight_layout()
#    plt.close()
#    fig.savefig(os.path.join(OUTER_Path_to_Folder,str(TimeSteps[j])+"_TimeOfPeakHist"+".pdf"), dpi=fig.dpi, bbox_inches='tight')

##BigUpSlope
#    histbins = numpy.linspace(min(BigUpSlope[j]),max(BigUpSlope[j]),20)
#    hist, bins = numpy.histogram(BigUpSlope[j], bins=histbins, density=True)
#    width = 0.7 * (bins[1] - bins[0])
#    center = (bins[:-1] + bins[1:]) / 2
#    
#    mean, STD = norm.fit(BigUpSlope[j])
#    x = numpy.linspace(min(BigUpSlope[j]), max(BigUpSlope[j]), 1000)
#    p = norm.pdf(x, mean, STD)
#    
#    plt.rcParams['figure.figsize'] = [10, 10]
#    fig, ax = plt.subplots()
#    ax.bar(center, hist, align='center', width=width)
#    plt.plot(x, p)
#    plt.tight_layout()
#    plt.close()
#    fig.savefig(os.path.join(OUTER_Path_to_Folder,str(TimeSteps[j])+"_UpSlopeHist"+".pdf"), dpi=fig.dpi, bbox_inches='tight')

##BigDownSlope
#    histbins = numpy.linspace(min(BigDownSlope[j]),max(BigDownSlope[j]),20)
#    hist, bins = numpy.histogram(BigDownSlope[j], bins=histbins, density=True)
#    width = 0.7 * (bins[1] - bins[0])
#    center = (bins[:-1] + bins[1:]) / 2
#    
#    mean, STD = norm.fit(BigDownSlope[j])
#    x = numpy.linspace(min(BigDownSlope[j]), max(BigDownSlope[j]), 1000)
#    p = norm.pdf(x, mean, STD)
#    
#    plt.rcParams['figure.figsize'] = [10, 10]
#    fig, ax = plt.subplots()
#    ax.bar(center, hist, align='center', width=width)
#    plt.plot(x, p)
#    plt.tight_layout()
#    plt.close()
#    fig.savefig(os.path.join(OUTER_Path_to_Folder,str(TimeSteps[j])+"_DownSlopeHist"+".pdf"), dpi=fig.dpi, bbox_inches='tight')

##BigEndTime
#    histbins = numpy.linspace(min(BigEndTime[j]),max(BigEndTime[j]),20)
#    hist, bins = numpy.histogram(BigEndTime[j], bins=histbins, density=True)
#    width = 0.7 * (bins[1] - bins[0])
#    center = (bins[:-1] + bins[1:]) / 2
#    
#    mean, STD = norm.fit(BigEndTime[j])
#    x = numpy.linspace(min(BigEndTime[j]), max(BigEndTime[j]), 1000)
#    p = norm.pdf(x, mean, STD)
#    
#    plt.rcParams['figure.figsize'] = [10, 10]
#    fig, ax = plt.subplots()
#    ax.bar(center, hist, align='center', width=width)
#    plt.plot(x, p)
#    plt.tight_layout()
#    plt.close()
#    fig.savefig(os.path.join(OUTER_Path_to_Folder,str(TimeSteps[j])+"_DownSlopeHist"+".pdf"), dpi=fig.dpi, bbox_inches='tight')

##BigAUC
#    histbins = numpy.linspace(min(BigAUC[j]),max(BigAUC[j]),20)
#    hist, bins = numpy.histogram(BigAUC[j], bins=histbins, density=True)
#    width = 0.7 * (bins[1] - bins[0])
#    center = (bins[:-1] + bins[1:]) / 2
#    
#    mean, STD = norm.fit(BigAUC[j])
#    x = numpy.linspace(min(BigAUC[j]), max(BigAUC[j]), 1000)
#    p = norm.pdf(x, mean, STD)
#    
#    plt.rcParams['figure.figsize'] = [10, 10]
#    fig, ax = plt.subplots()
#    ax.bar(center, hist, align='center', width=width)
#    plt.plot(x, p)
#    plt.tight_layout()
#    plt.close()
#    fig.savefig(os.path.join(OUTER_Path_to_Folder,str(TimeSteps[j])+"_DownSlopeHist"+".pdf"), dpi=fig.dpi, bbox_inches='tight')

#################################################################################
#Makes Aspect Graphs
#################################################################################
    size = 32
    capsizes = [50, 21, 50, 21, 50, 21, 50, 21, 50, 21]
    axessize = 2

    plt.rcParams['xtick.labelsize'] = size
    plt.rcParams['ytick.labelsize'] = size
    plt.rcParams['axes.labelsize'] = size
    plt.rcParams['figure.figsize'] = [12 , 8]
    plt.rcParams['xtick.major.size'] = axessize*4
    plt.rcParams['xtick.major.width'] = axessize
    plt.rcParams['xtick.minor.width'] = axessize
    plt.rcParams['ytick.major.width'] = axessize
    plt.rcParams['ytick.minor.width'] = axessize
    plt.rc('axes', linewidth=axessize)
    
    fig, ax = plt.subplots()
    ax.yaxis.set_major_locator(ticker.MaxNLocator(6))
    holder = []
    holder1 = []
    for i in range(len(TimeSteps)):
        holder.append(numpy.mean(BigPeakVirus[i][:]))
        holder1.append(numpy.std(BigPeakVirus[i][:]))
    plt.errorbar(TimeSteps.astype(numpy.float), holder, yerr=holder1, fmt = "-o", ecolor="black", capsize=10, linewidth = 3)
    plt.xlabel("TimeSteps")
    plt.ylabel("Peak Virus")
    plt.xscale("log")
    ax.set_xticks(TimeSteps.astype(numpy.float))
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    plt.xticks(rotation=-20)
    plt.tight_layout()
    plt.close()
    fig.savefig(os.path.join(OUTEROUTER_Path_to_Folder, "PeakViralTiter"+".pdf"), dpi=fig.dpi, bbox_inches='tight')
    
    fig, ax = plt.subplots()
    ax.yaxis.set_major_locator(ticker.MaxNLocator(6))
    holder = []
    holder1 = []
    for i in range(len(TimeSteps)):
        holder.append(numpy.mean(BigTimeOfPeak[i][:])/24)
        holder1.append(numpy.std(BigTimeOfPeak[i][:])/24)
    plt.errorbar(TimeSteps.astype(numpy.float), holder, yerr=holder1, fmt = "-o", ecolor="black", capsize=10, linewidth = 3)
    plt.xlabel("TimeSteps")
    plt.ylabel("Time of Peak (day)")
    plt.xscale("log")
    ax.set_xticks(TimeSteps.astype(numpy.float))
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    plt.xticks(rotation=-20)
    plt.tight_layout()
    plt.close()
    fig.savefig(os.path.join(OUTEROUTER_Path_to_Folder, "TimeofPeakViralTiter"+".pdf"), dpi=fig.dpi, bbox_inches='tight')
    
    fig, ax = plt.subplots()
    ax.yaxis.set_major_locator(ticker.MaxNLocator(6))
    holder = []
    holder1 = []
    for i in range(len(TimeSteps)):
        holder.append(numpy.mean(BigUpSlope[i][:])*24)
        holder1.append(numpy.std(BigUpSlope[i][:])*24)
    plt.errorbar(TimeSteps.astype(numpy.float), holder, yerr=holder1, fmt = "-o", ecolor="black", capsize=10, linewidth = 3)
    plt.xlabel("TimeSteps")
    plt.ylabel("Growth Rate (virus/day)")
    plt.xscale("log")
    ax.set_xticks(TimeSteps.astype(numpy.float))
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    plt.xticks(rotation=-20)
    plt.tight_layout()
    plt.close()
    fig.savefig(os.path.join(OUTEROUTER_Path_to_Folder, "UpSlope"+".pdf"), dpi=fig.dpi, bbox_inches='tight')
    
    fig, ax = plt.subplots()
    ax.yaxis.set_major_locator(ticker.MaxNLocator(6))
    holder = []
    holder1 = []
    holder2 = []
    holder3 = []
    for i in range(len(TimeSteps)):
        holder.append(numpy.mean(BigDownSlope[i][:])*24)
        holder1.append(numpy.std(BigDownSlope[i][:])*24)
        holder2.append(numpy.round(numpy.mean(BigDownSlope[i][:]),3)*24)
        holder3.append(numpy.round(numpy.std(BigDownSlope[i][:]),3)*24)
    plt.errorbar(TimeSteps.astype(numpy.float), holder, yerr=holder1, fmt = "-o", ecolor="black", capsize=10, linewidth = 3)
    plt.xlabel("TimeSteps")
    plt.ylabel("Decay Rate (virus/day)")
    plt.xscale("log")
    ax.set_xticks(TimeSteps.astype(numpy.float))
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    plt.xticks(rotation=-20)
    plt.tight_layout()
    plt.close()
    fig.savefig(os.path.join(OUTEROUTER_Path_to_Folder, "DownSlope"+".pdf"), dpi=fig.dpi, bbox_inches='tight')
#    fig, ax = plt.subplots()
#    ax.yaxis.set_major_locator(ticker.MaxNLocator(6))
#    plt.errorbar(TimeSteps.astype(numpy.float), holder2, yerr=holder3, fmt = "-o", ecolor="black", capsize=10, linewidth = 3)
#    plt.xlabel("TimeSteps")
#    plt.ylabel("Decay Rate (virus/day)")
#    plt.xscale("log")
#    ax.set_xticks(TimeSteps.astype(numpy.float))
#    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
#    plt.xticks(rotation=-20)
#    plt.tight_layout()
#    plt.close()
#    fig.savefig(os.path.join(OUTEROUTER_Path_to_Folder, "DownSlopeRound"+".pdf"), dpi=fig.dpi, bbox_inches='tight')

    fig, ax = plt.subplots()
    ax.yaxis.set_major_locator(ticker.MaxNLocator(6))
    holder = []
    holder1 = []
    for i in range(len(TimeSteps)):
        holder.append(numpy.mean(BigEndTime[i][:])/24)
        holder1.append(numpy.std(BigEndTime[i][:])/24)
    plt.errorbar(TimeSteps.astype(numpy.float), holder, yerr=holder1, fmt = "-o", ecolor="black", capsize=10, linewidth = 3)
    plt.xlabel("TimeSteps")
    plt.ylabel("Infection Duration (day)")
    plt.xscale("log")
    ax.set_xticks(TimeSteps.astype(numpy.float))
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    plt.xticks(rotation=-20)
    plt.tight_layout()
    plt.close()
    fig.savefig(os.path.join(OUTEROUTER_Path_to_Folder, "InfectionDuration"+".pdf"), dpi=fig.dpi, bbox_inches='tight')

    fig, ax = plt.subplots()
    ax.yaxis.set_major_locator(ticker.MaxNLocator(6))
    holder = []
    holder1 = []
    for i in range(len(TimeSteps)):
        holder.append(numpy.mean(BigAUC[i][:])/24)
        holder1.append(numpy.std(BigAUC[i][:])/24)
    plt.errorbar(TimeSteps.astype(numpy.float), holder, yerr=holder1, fmt = "-o", ecolor="black", capsize=10, linewidth = 3)   
    plt.xlabel("TimeSteps")
    plt.ylabel("AUC (log(virus)*day)")
    plt.xscale("log")
    ax.set_xticks(TimeSteps.astype(numpy.float))
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    plt.xticks(rotation=-20)
    plt.tight_layout()
    plt.close()
    fig.savefig(os.path.join(OUTEROUTER_Path_to_Folder, "AUC"+".pdf"), dpi=fig.dpi, bbox_inches='tight')

#################################################################################
    with open(os.path.join(OUTEROUTER_Path_to_Folder,"PeakVirus.txt"), 'w') as outfile:
        for i in range(len(TimeSteps)):
            print(str(numpy.mean(BigPeakVirus[i][:]))+","+str(numpy.std(BigPeakVirus[i][:])),file=outfile)
        
    with open(os.path.join(OUTEROUTER_Path_to_Folder,"TimeOfPeakVirus.txt"), 'w') as outfile:
        for i in range(len(TimeSteps)):
            print(str(numpy.mean(BigTimeOfPeak[i][:]))+","+str(numpy.std(BigTimeOfPeak[i][:])),file=outfile)
            
    with open(os.path.join(OUTEROUTER_Path_to_Folder,"UpSlope.txt"), 'w') as outfile:
        for i in range(len(TimeSteps)):
            print(str(numpy.mean(BigUpSlope[i][:]))+","+str(numpy.std(BigUpSlope[i][:])),file=outfile)
    with open(os.path.join(OUTEROUTER_Path_to_Folder,"DownSlope.txt"), 'w') as outfile:
        for i in range(len(TimeSteps)):
            print(str(numpy.mean(BigDownSlope[i][:]))+","+str(numpy.std(BigDownSlope[i][:])),file=outfile)
    with open(os.path.join(OUTEROUTER_Path_to_Folder,"InfectionDuration.txt"), 'w') as outfile:
        for i in range(len(TimeSteps)):
            print(str(numpy.mean(BigEndTime[i][:]))+","+str(numpy.std(BigEndTime[i][:])),file=outfile)
    with open(os.path.join(OUTEROUTER_Path_to_Folder,"AUC.txt"), 'w') as outfile:
        for i in range(len(TimeSteps)):
            print(str(numpy.mean(BigAUC[i][:]))+","+str(numpy.std(BigAUC[i][:])),file=outfile)

##################################################################################
##Test if the distribution of the measurable is normal
##################################################################################
#    BigArray = [BigPeakVirus, BigTimeOfPeak, BigUpSlope, BigDownSlope, BigEndTime, BigAUC]
#    BigStrArray = ["BigPeakVirus", "BigTimeOfPeak", "BigUpSlope", "BigDownSlope", "BigEndTime", "BigAUC"]
#    
#    for q in range(len(BigArray)):  
#        #    Shapiro-Wilk
#        SStat, SP = shapiro(BigArray[q][j])
#        #    D'Agostino's
#        DStat, DP = normaltest(BigArray[q][j])
#        #    Anderson-Darling
#        result = anderson(BigArray[q][j])
#        
#        with open(os.path.join(OUTER_Path_to_Folder,str(TimeSteps[j])+"_Normality"+".txt"),'a') as outfile:
#        
#            print(str(BigStrArray[q]),file=outfile)
#        
#        #    Shapiro-Wilk
#            print("Shapiro-Wilk",file=outfile)
#            print("Statistics= "+"{0:.3g}".format(SStat)+", p= "+"{0:.3g}".format(SP),file=outfile)
#            SAlpha = 0.05
#            if SP > SAlpha:
#	            print('Sample looks Gaussian (fail to reject H0)',file=outfile)
#            else:
#	            print('Sample does not look Gaussian (reject H0)',file=outfile)
#	            
#        #    D'Agostino's
#            print("D'Agostino's",file=outfile)
#            print("Statistics= "+"{0:.3g}".format(DStat)+", p= "+"{0:.3g}".format(DP),file=outfile)
#            DAlpha = 0.05
#            if DP > DAlpha:
#	            print('Sample looks Gaussian (fail to reject H0)',file=outfile)
#            else:
#	            print('Sample does not look Gaussian (reject H0)',file=outfile)
#	            
#        #    Anderson-Darling
#            print("Anderson-Darling",file=outfile)
#            print("Statistic: "+"{0:.3g}".format(result.statistic),file=outfile)
#            p = 0
#            for i in range(len(result.critical_values)):
#	            sl, cv = result.significance_level[i], result.critical_values[i]
#	            if result.statistic < result.critical_values[i]:
#		            print("{0:.3g}".format(sl)+":"+"{0:.3g}".format(cv)+", data looks normal (fail to reject H0)",file=outfile)
#	            else:
#		            print("{0:.3g}".format(sl)+":"+"{0:.3g}".format(cv)+", data does not looks normal (reject H0)",file=outfile)
#            
#            print("\n",file=outfile)

#################################################################################
#Plots all the runs and the median curves
#################################################################################
directory = os.path.join("".join(OUTERPATH+ "_" +"AllOnOne"))
if not os.path.exists(directory):
    os.makedirs(directory)
AllOnOne_Path_to_Folder = os.path.abspath(directory)
#All on one
size = 32
axessize = 2
days = 11
plt.rcParams['xtick.labelsize'] = size
plt.rcParams['ytick.labelsize'] = size
plt.rcParams['axes.labelsize'] = size
plt.rc('legend', fontsize=size*(6/10))
plt.rc('legend', title_fontsize=size*(6/10))
plt.rcParams['figure.figsize'] = [10, 10]
plt.rcParams['xtick.major.size'] = axessize*4
plt.rcParams['xtick.major.width'] = axessize
plt.rcParams['xtick.minor.width'] = axessize
plt.rcParams['ytick.major.width'] = axessize
plt.rcParams['ytick.minor.width'] = axessize
plt.rc('axes', linewidth=axessize)
fig, ax = plt.subplots()
plt.yscale("log")
plt.ylim(numpy.power(10,1), numpy.power(10,15))
plt.xlim(0, days*24)
#colors = ['#762a83', '#af8dc3', '#e7d4e8', '#f7f7f7', '#d9f0d3', '#7fbf7b', '#1b7837']
colors = ['#762a83', '#af8dc3', '#e7d4e8', '#000000', '#d9f0d3', '#7fbf7b', '#1b7837']
#alphas = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.05, 0.05]
alphas = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
#alphas = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
for j in range(len(TimeSteps)):
#    if(j%2==1):
    for i in ReallyBigVirusArray[j]:
        plt.plot(i, color = colors[j], alpha = alphas[j])
    plt.plot(BigMedianVirusArray[j], color = colors[j], label = TimeSteps[j])
#plt.hlines(numpy.power(10,7), 0, days*24, colors='k', linestyles='dashed')
plt.xticks(numpy.linspace(0, days*24, days), numpy.linspace(0, days*24, days, dtype = int)//24)
plt.xticks(rotation=-20)
plt.tight_layout()
plt.ylabel("Median Virus Titer")
plt.xlabel("Time (day)")
plt.legend(title="time step")
plt.close()
fig.savefig(os.path.join(AllOnOne_Path_to_Folder,"AllOnOne_Median_VirusVsTime"+".pdf"), dpi=fig.dpi, bbox_inches='tight')

###############################################################################
print("Data Analysis Complete")
