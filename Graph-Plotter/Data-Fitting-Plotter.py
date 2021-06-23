import numpy as np
import matplotlib.pyplot as plt
import os
import gc
import pathlib
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from scipy.stats import shapiro

SCTime = [2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 15, 18]
SC = [ 0.00007109709432312438, 0.0002586416205275976, 0.0019073187599725216, 0.039137456019803916, 0.048733846139805366, 0.16885472875899837, 0.2618119214401843, 0.31815994481219106, 0.3961718357595989, 0.2493592004984161, 0.4262158829015338, 1.2152232910634688]

MCTime = [12, 18, 24, 30, 36, 48, 60, 72, 84, 96, 108, 120]
MC = np.array([819.4679979024306, 7759.590246392491, 2905620.46821282, 73765857.98521078, 999437728.9462978, 6640163884.645269, 6709615082.419942, 5663087196.46287, 1040447862.7422994, 1051330177.0262353, 54967688.493010186, 1403746.5929611742])

#MCError = np.empty([len(MC)])
#for i in range(len(MC)):   
#    power = np.around(MC[i],-1*int(np.floor(np.log10(MC[i]))+1))
#    if power == 0.0:
#        MCError[i] = np.power(10.0,int(np.floor(np.log10(MC[i]))))
#    else:
#        MCError[i] = np.around(MC[i],-1*int(np.floor(np.log10(MC[i]))+1))
#print(MC)
#print(MCError)
#MCError = MCError*0.5
#print(MCError)

MCError = MC*0.5

#SSRArray = np.empty([625],dtype=tuple )

#with open("Results_SSR.txt") as infile:
#    index = 0
#    for line in infile:
#        Data = line
#        Data = Data.replace(" ","")
#        Data = Data.replace(":",",")
#        Data = Data.replace("\n","")
#        Data = Data.split(",")
#        del Data[0]
#        del Data[1]
#        del Data[2]
#        del Data[3]
#        del Data[4]
#        SSRArray[index] = np.asarray(Data, dtype=np.float)
#        index = index + 1

#Max = 0.0
#for i in range(625):
#    if(SSRArray[i][0] > Max):
#        Max = SSRArray[i][0]

##betaarray = [1.0, 1.5, 2.0, 2.5, 3.0]
##for j in range(len(betaarray)):
##    fig = plt.figure()
##    ax = fig.add_subplot(111, projection='3d')
##    for i in range(625):
##        if(SSRArray[i][1] == betaarray[j]):
##            cnum = np.log10(float(SSRArray[i][0]))/np.log10(Max)
##            ax.scatter(SSRArray[i][3], SSRArray[i][4], SSRArray[i][2], marker='o', alpha=(1-cnum))
##    ax.set_xlabel('TauI')
##    ax.set_ylabel('TauE')
##    ax.set_zlabel('Rho')
##    ax.set_title(str(betaarray[j]))
##    plt.show()

#x = np.empty([5,5])
#y = np.empty([5,5])
#z = np.empty([5,5])

#u = np.empty([25])
#v = np.empty([25])
#w = np.empty([25])

#idx = 0
#idy = 0
#index = 0
#for i in range(625):
#    if((SSRArray[i][3] == 6.0) & (SSRArray[i][4] == 9.0)):
#        x[idx][idy] = SSRArray[i][2]
#        y[idx][idy] = SSRArray[i][1]
#        z[idx][idy] = SSRArray[i][0]

#        u[index] = SSRArray[i][2]
#        v[index] = SSRArray[i][1]
#        w[index] = SSRArray[i][0]
#        
#        index = index + 1
#        idx = idx + 1
#        if(idx == 5):
#            idy = idy + 1
#            idx = 0
#        if(idy == 5):
#            idy = 0

#minid = np.where(w==1.5552109355255736e+20)[0]

#plt.rcParams['figure.figsize'] = [10, 10]
#fig = plt.figure()
#ax = fig.gca(projection='3d')

##ax.plot_surface(x, y, z, alpha=0.8, cmap=cm.coolwarm)
##ax.plot_trisurf(u, v, w, alpha=0.5, color="green")
#ax.plot(u, v, w, "o", alpha=1.0, color="black")
#ax.plot(u[minid], v[minid], w[minid], "o", alpha=1.0, color="magenta")
##ax.scatter(u, v, w, c="k", cmap=cm.coolwarm, alpha=1.0)

##ax.contourf(x, y, z, zdir='x', offset=x.max(), cmap=cm.coolwarm, alpha=0.5)
##ax.contourf(x, y, z, zdir='y', offset=y.max(), cmap=cm.coolwarm, alpha=0.5)
##ax.contourf(x, y, z, zdir='z', offset=z.min(), cmap=cm.coolwarm, alpha=0.5)

#ax.set_xlabel(r'$\beta$')
#ax.set_ylabel(r'$\rho$')
#ax.set_zlabel('SSR')

#ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
#ax.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
#ax.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))

#ax.set_xlim(x.min(), x.max())
#ax.set_ylim(y.min(), y.max())
#ax.set_zlim(z.min(), z.max())

#ax.view_init(30, -135)

#plt.axis('off')

##plt.show()
#fig.savefig("4.pdf", dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)

#exit(0)
##############################################################
#Initial  
#beta = 54.21
#rho = 3000.0
#TauI = 26.24
#TauE = 15.16
#c = 0.25



beta = 53.64218701755928
rho = 3005.675224284383
TauI = 26.188840368473215
TauE = 15.874079384166595
c = 0.2515931715378636

print("Run:"+str(beta)+","+str(rho)+","+str(TauI)+","+str(TauE)+","+str(c))

NumberOfRuns = 10
VirusData = np.empty([NumberOfRuns],dtype=tuple )
MedianData = []
MedianCD = []
LongestTime = 0

#with open("/home/baylorfain/Documents/CompSci-Convergance-Testing/Data-Fitting/Run:"+str(beta)+","+str(rho)+","+str(TauI)+","+str(TauE)+","+str(c)+".txt") as infile:
#with open("/home/baylorfain/Documents/Data-Fitting-Results/SSR/10^-4/Run:"+str(beta)+","+str(rho)+","+str(TauI)+","+str(TauE)+","+str(c)+".txt") as infile:
with open("Run:"+str(beta)+","+str(rho)+","+str(TauI)+","+str(TauE)+","+str(c)+".txt") as infile:
    index = 0
    for line in infile:
        Data = line
        Data = Data.split(",")
        del Data[-1]
        VirusData[index] = np.asarray(Data, dtype=np.float)
        if(LongestTime < len(VirusData[index])):
            LongestTime = len(VirusData[index])
        index = index + 1
        if index == NumberOfRuns:
            break

    for j in range(LongestTime):
        InbetweenArray = []
        for i in range(NumberOfRuns):
            try:
                InbetweenArray.append(VirusData[i][j])
            except IndexError:
                InbetweenArray.append(0.0) #maybr toggle this
    
        MedianData.append(np.median(InbetweenArray))

        for time in MCTime:
            if j == time:
                SStat, SP = shapiro(InbetweenArray)
                SAlpha = 0.05
                if SP > SAlpha:
                    print('Normal')
                else:
                    print('Not Normal')

size = 32
axessize = 2
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

ax.plot(np.log10(MedianData), color="blue", label = "Simulation")
for i in range(NumberOfRuns):
    ax.plot(np.log10(VirusData[i]), color="blue", alpha=0.1)
ax.errorbar(MCTime, np.log10(MC), yerr=0.5, fmt="o", color="green", label = "Experimental")

xmin = 6
xmax = MCTime[-1]+6
ymin = 2
ymax = 10
#If x is longer than y
ax.set_aspect((xmax-xmin)/(ymax-ymin), adjustable='box')
#If y is longer than x
#ax.set_aspect((ymax-ymin)(xmax-xmin), adjustable='box')
plt.xlim(xmin,xmax)
plt.ylim(ymin,ymax)
plt.xticks(MCTime, MCTime)
plt.xticks(rotation=-45)
ticks = np.arange(ymin, ymax+1, 1)
ticklabels = [r"$10^{{{}}}$".format(tick) for tick in ticks]
plt.yticks(ticks, ticklabels)
plt.tight_layout()
plt.xlabel('Time (h)')
plt.ylabel('Virus (PFU/ml or RNA/ml)')
plt.legend(title="Data Source")
plt.close()
fig.savefig("DataFit"+".pdf", dpi=fig.dpi, bbox_inches='tight')
#plt.show()

###################################################

#MOI = [-3.0]

#for i in range(len(MOI)):
#    CELL2CELLPATH = r"/media/baylor/My Passport/BaylorFain/500Runs/Cell2Cell/607_0-CELL2CELL_"+str(MOI[i])+"-MOI"
#    CELLFREEPATH = r"/media/baylor/My Passport/BaylorFain/500Runs/CellFree/607_0-FREECELL_"+str(MOI[i])+"-MOI"
#    SAVEPATH = r"/home/baylor/Documents/Research/Presentations/All_Presentations/SIAM_TX-LA/2020"

#    cell2cellvirus = []
#    cellfreevirus = []
#    with open(os.path.join(CELL2CELLPATH,"PerTimeStep.txt")) as cell2cellinfile:
#        with open(os.path.join(CELLFREEPATH,"PerTimeStep.txt")) as cellfreeinfile:
#            for line in cell2cellinfile:
#                line = line.split(",")
#                cell2cellvirus.append(float(line[5]))
#                
#            for line in cellfreeinfile:
#                line = line.split(",")
#                cellfreevirus.append(float(line[5]))

#    size = 32
#    plt.rcParams['xtick.labelsize'] = size
#    plt.rcParams['ytick.labelsize'] = size
#    plt.rcParams['axes.labelsize'] = size
#    plt.rcParams['figure.figsize'] = [10, 10]
#    fig, ax = plt.subplots()
#    plt.yscale("log")
#    DaysOut = 24*40
#    Cell2CellTime = (np.linspace(0,len(cell2cellvirus[0:DaysOut]),len(cell2cellvirus[0:DaysOut])))/24
#    CellFreeTime = (np.linspace(0,len(cellfreevirus[0:DaysOut]),len(cellfreevirus[0:DaysOut])))/24
##    plt.plot(Cell2CellTime[0:DaysOut], cell2cellvirus[0:DaysOut],"-", label = "Cell to Cell")
#    plt.plot(CellFreeTime[0:DaysOut], cellfreevirus[0:DaysOut],"-", label = "Cell Free")
#    plt.xlabel("Time (days)")
#    plt.ylabel("Amount of Virus")
#    plt.legend(prop={'size': size*(7/10)})
#    plt.tight_layout()
#    plt.close()
#    fig.savefig(os.path.join(SAVEPATH,"MOI_{}_".format(MOI[i]) + "CellToCell_CellFree_ViralTitter"+".pdf"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)



