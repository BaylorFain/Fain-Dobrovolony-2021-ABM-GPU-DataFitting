#import pylab as plt
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy
from scipy import interpolate
from scipy import ndimage
import os
import gc
from mpl_toolkits import mplot3d
import seaborn as sns
import pandas as pd


string = r"/media/baylor/TOSHIBA EXT/Comp-Sci-Data/Newstuff/607_0-CELLFREE_-4.0-MOI_0.005000-timestep"

NumberOfLayers = 607

Nx=(2*NumberOfLayers)-1
Ny=(2*NumberOfLayers)-1

s = (2.0/3.0)
RadiusScale = 0
for i in range(NumberOfLayers):
    if i == 0:
        RadiusScale = RadiusScale + 1
    else:
        if (i)%2 == 1:
            RadiusScale = RadiusScale + 1
        else:
            RadiusScale = RadiusScale + 2
RadiusOfCircle = s*RadiusScale

count = 0
for i in range(NumberOfLayers):
    count = count + i
N=(count)*6+1

coord = numpy.empty([N,3],dtype=float )
percyclecoord = numpy.empty([N,3],dtype=tuple)

percyclecoord.fill(0)
coord.fill(0)

percyclecoord[0] = [0,0,0]
for j in range(NumberOfLayers):
    for i in range(2*j):
        if i < j:
            temp = i
        percyclecoord[i+(j-1)*j+1][0] =  -temp-1
        percyclecoord[i+(j-1)*j+1][1] =   temp+j-i #Which is just -x-z
        percyclecoord[i+(j-1)*j+1][2] =  -j+1+i

c0 = [percyclecoord[0][0], percyclecoord[0][1], percyclecoord[0][2]]
coord[0][2] = c0[2]
coord[0][1] = c0[1]
coord[0][0] = c0[0]

count = 0
for j in range(int(N/3)):
    for i in range(3):
        coord[(i+0)%3+3*j+1][2] = percyclecoord[j+1][i]+c0[i]
        coord[(i+1)%3+3*j+1][1] = percyclecoord[j+1][i]+c0[i]
        coord[(i+2)%3+3*j+1][0] = percyclecoord[j+1][i]+c0[i]

hi = coord[0][0]
vi = coord[0][2]
xmin = numpy.Inf
for i in range(len(coord)):
    xcoord = coord[i][0]
    if coord[i][0] < xmin:
        xmin = coord[i][0]
    ycoord = (2.0 * numpy.sin(numpy.radians(60)) * (coord[i][1] - coord[i][2]) /3.0)+vi
    dist = numpy.sqrt((xcoord-hi)**2+(ycoord-vi)**2)
    if dist >= RadiusOfCircle:
        coord[i][0] = 5000.0
        coord[i][1] = 0.0
        coord[i][2] = 0.0

QRIndexing = numpy.empty([(2*NumberOfLayers)-1,(2*NumberOfLayers)-1],dtype=float)
QRIndexing.fill(0.0)
for i in range(len(coord)):
    if coord[i][0] != 5000:
        QRIndexing[int(coord[i][2])-int(xmin)][int(coord[i][0])-int(xmin)] = 1.0
Index = numpy.where(QRIndexing != 0)
CoordanateNested = numpy.empty([(2*NumberOfLayers)-1,(2*NumberOfLayers)-1],dtype=tuple)
CoordanateNested.fill(0)
for i in range(len(Index[0])):
    CoordanateNested[Index[0][i]][Index[1][i]] = [Index[0][i]+xmin,-(Index[1][i]+xmin+Index[0][i]+xmin),Index[1][i]+xmin]

QRIndexing1 = numpy.empty([(2*NumberOfLayers)-1,(2*NumberOfLayers)-1],dtype=float)
QRIndexing1.fill(1)
Index1 = numpy.where(QRIndexing1 != 0)
CoordanateNested1 = numpy.empty([(2*NumberOfLayers)-1,(2*NumberOfLayers)-1],dtype=tuple)
CoordanateNested1.fill(0)
for i in range(len(Index1[0])):
    CoordanateNested1[Index1[0][i]][Index1[1][i]] = [1.0*Index1[0][i]+1.0*xmin,-(1.0*Index1[1][i]+1.0*xmin+1.0*Index1[0][i]+1.0*xmin),1.0*Index1[1][i]+1.0*xmin]


INITIALPATH = os.path.join(string).replace('\\','//')

directory = os.path.join("".join(INITIALPATH + "/VirusPhotos"))
if not os.path.exists(directory):
    os.makedirs(directory)
Path_to_Folder = os.path.abspath(directory)

TimeArray = []
HealthyArray = []
EclipseArray = []
InfectedArray = []
DeadArray = []
VirusArray = []

print("Biginning Sorting")
with open(os.path.join(INITIALPATH,"PerTimeStep.txt")) as outfile:
    for cellline in outfile:
        CellData = cellline
        CellData = CellData.split(",")
        del CellData[-1]
            
#        TimeArray.append(int(CellData[0]))
#        HealthyArray.append(int(CellData[1]))
#        EclipseArray.append(int(CellData[2]))
#        InfectedArray.append(int(CellData[3]))
#        DeadArray.append(int(CellData[4]))
        VirusArray.append(float(CellData[5]))
IndexOfPeakTime = numpy.argmax(VirusArray)

print("Biginning Max/Min Find")
with open(os.path.join(INITIALPATH,"virus_over_time.txt")) as VirusInfile:
    maxvirusdata = 0.0 
    minvirusdata = numpy.inf
    Time = 0
    for virusline in VirusInfile:
        VirusData = virusline
        VirusData = VirusData.split(",")
        del VirusData[Nx]
        for q in range(len(VirusData)):
            number = float(VirusData[q])
            if number>maxvirusdata:
                maxvirusdata = number
            elif (number<minvirusdata)&(number!=0.0):
                minvirusdata = number
        if (Time%Nx==0)&(VirusArray[Time//Nx]<=(0.7*VirusArray[IndexOfPeakTime]))&((Time//Nx)>IndexOfPeakTime):
            break
        Time = Time +1
levels = numpy.linspace(numpy.log10(minvirusdata), numpy.log10(maxvirusdata),50)

print("Printing Biginnings")
rgba_colors = numpy.zeros([Nx*Ny,4])
rgba_colors[:,0] = 0.5
rgba_colors[:,1] = 0.0
rgba_colors[:,2] = 0.5

hc = numpy.empty([Nx*Ny],dtype=float )
vc = numpy.empty([Nx*Ny],dtype=float )
hcv = numpy.empty([Nx, Ny],dtype=float )
vcv = numpy.empty([Nx, Ny],dtype=float )
VirusDataArray = numpy.empty([Nx, Ny],dtype=float )

hc.fill(0.0)
vc.fill(0.0)
hcv.fill(0.0)
vcv.fill(0.0)

with open(os.path.join(INITIALPATH,"virus_over_time.txt")) as VirusInfile:
    index = 0
    count = 0
    for virusline in VirusInfile:
        VirusData = virusline
        VirusData = VirusData.split(",")
        del VirusData[Nx]
        for i in range(len(VirusData)):
            if float(VirusData[i]) == 0.0:
                VirusDataArray[index][i] = numpy.nan
            else:
                VirusDataArray[index][i] = numpy.log10(float(VirusData[i]))
        index = index + 1
        
        if index == (Nx):
            for n in range(Nx):
                for i in range(Nx):
                    number = VirusDataArray[i][n]/numpy.log10(maxvirusdata)
                    if number < 0.0:
                        number = 0.0
                    rgba_colors[n+Nx*i,3] = number
                    if CoordanateNested[i][n] != 0:
                        hc[n+Nx*i] = s*(CoordanateNested[i][n][0] - 0.5*(CoordanateNested[i][n][1] + CoordanateNested[i][n][2]))
                        vc[n+Nx*i] = s*(numpy.sqrt(3)/2*(CoordanateNested[i][n][1] - CoordanateNested[i][n][2]))
                        
                    hcv[i][n] = s*(CoordanateNested1[i][n][0] - 0.5*(CoordanateNested1[i][n][1] + CoordanateNested1[i][n][2]))
                    vcv[i][n] = s*(numpy.sqrt(3)/2*(CoordanateNested1[i][n][1] - CoordanateNested1[i][n][2]))

            #Making pictures
            plt.rcParams['figure.figsize'] = [10, 10]
            fig, ax = plt.subplots()
            plt.axis("off")
            plt.axis('equal')
            if count == 0:
                ax.plot(hc,vc, color=(1.0,1.0,1.0), marker = ",", linestyle = "none", markersize = 0.1)
                [ymin, ymax] = ax.get_ylim()
                [xmin, xmax] = ax.get_xlim()
            #Full Dish
            ax.set_ylim(ymin, ymax)
            ax.set_xlim(xmin, xmax)
            #Zoomed In
#            zoomx = xmin+(442/769)*2*xmax
#            zoomy = ymin+((770-363)/770)*2*ymax
#            ax.set_ylim(zoomy-50, zoomy+50)
#            ax.set_xlim(zoomx-50, zoomx+50)



            ax.contourf(hcv,vcv,VirusDataArray, levels=levels)
            #, cmap=plt.get_cmap("Purples")



            plt.tight_layout()
            plt.close()
            fig.savefig(os.path.join(Path_to_Folder,"virus"+"{}".format(count)+".png"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)

            plt.cla()
            plt.clf()
            plt.close('all')
            gc.collect()

            index = 0
            print(count)
            count = count + 1




