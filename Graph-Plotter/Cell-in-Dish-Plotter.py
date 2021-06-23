import numpy
import datetime
import matplotlib.pyplot as plt
import os
import gc

os.system('cls' if os.name=='nt' else 'clear')



NumberOfLayers = 607



#string = r"/media/baylor/TOSHIBA EXT/Comp-Sci-Data/Plaque-Movie/Save-every-hour/607_0-CELLFREE_-4.0-MOI_0.005000-timestep"
string = r"/media/baylor/TOSHIBA EXT/Comp-Sci-Data/Newstuff/607_0-CELLFREE_-4.0-MOI_0.005000-timestep"



INITIALPATH = os.path.join(string).replace('\\','//')
directory = os.path.join("".join(INITIALPATH + "/CellPhotos"))
if not os.path.exists(directory):
    os.makedirs(directory)
Path_to_Folder = os.path.abspath(directory) # Figures out the absolute path for you in case your working directory moves around.

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
    CoordanateNested[Index[0][i]][Index[1][i]] = [1.0*Index[0][i]+1.0*xmin,-(1.0*Index[1][i]+1.0*xmin+1.0*Index[0][i]+1.0*xmin),1.0*Index[1][i]+1.0*xmin]



Nx=(2*NumberOfLayers)-1
Ny=(2*NumberOfLayers)-1

CellDataArray = numpy.empty([Nx, Ny],dtype=tuple )
VirusDataArray = numpy.empty([Nx, Ny],dtype=tuple ) 

hx = []
hy = []
ex = []
ey = []
ix = []
iy = []
dx = []
dy = []

#    maxdata1 = 0.0 
#    mindata1 = numpy.inf
#    for i in range(N*NumberOfPics):
#        number = float(VirusData[i])
#        if number>maxdata1:
#            maxdata1 = number
#        elif number<mindata1:
#            mindata1 = number

print("Printing Biginnings")

with open(os.path.join(INITIALPATH,"cells_over_time.txt")) as CellInfile:
    index = 0    
    count = 0
    for cellline in CellInfile:
        CellData = cellline
        CellData = CellData.split(",")
        del CellData[Nx]
        for i in range(len(CellData)):
            CellDataArray[index][i] = CellData[i]
        index = index + 1

        if index == (Nx):
            for n in range(Nx):
                for i in range(Nx):
                    if CoordanateNested[i][n] != 0:
                        hc = s*(CoordanateNested[i][n][0] - 0.5*(CoordanateNested[i][n][1] + CoordanateNested[i][n][2]))
                        vc = s*(numpy.sqrt(3)/2*(CoordanateNested[i][n][1] - CoordanateNested[i][n][2]))

                    c = CellDataArray[i][n]
                    if c == "h":
                        hx.append(hc)
                        hy.append(vc)
                    elif c == "e":
                        ex.append(hc)
                        ey.append(vc)
                    elif c == "i":
                        ix.append(hc)
                        iy.append(vc)
                    elif c == "d":
                        dx.append(hc)
                        dy.append(vc)
                    elif c == "o":
                        "Nothing"
                        
            plt.rcParams['figure.figsize'] = [10, 10]
            fig, ax = plt.subplots()
            plt.axis("off")
            ax.set_aspect("equal")
            
            #Zoomed
            ax.plot(hx, hy, color = (0.0, 1.0, 0.0), alpha=0.5, marker = "H", linestyle = "none", markersize = 7, markeredgecolor = (0.0, 1.0, 0.0))
            ax.plot(ex, ey, color = (0.0, 1.0, 1.0), alpha=0.5, marker = "H", linestyle = "none", markersize = 7, markeredgecolor = (0.0, 1.0, 1.0))
            ax.plot(ix, iy, color = (1.0, 0.0, 0.0), alpha=0.5, marker = "H", linestyle = "none", markersize = 7, markeredgecolor = (1.0, 0.0, 0.0))
            ax.plot(dx, dy, color = (0.0, 0.0, 0.0), alpha=0.5, marker = "H", linestyle = "none", markersize = 7, markeredgecolor = (0.0, 0.0, 0.0))
            #Full Dish
#            ax.plot(hx, hy, color = (0.0, 1.0, 0.0), alpha=1.0, marker = ".", linestyle = "none", markersize = 1, markeredgewidth = 0.5)
#            ax.plot(ex, ey, color = (0.0, 1.0, 1.0), alpha=1.0, marker = ".", linestyle = "none", markersize = 1, markeredgewidth = 0.5)
#            ax.plot(ix, iy, color = (1.0, 0.0, 0.0), alpha=1.0, marker = ".", linestyle = "none", markersize = 1, markeredgewidth = 0.5)
#            ax.plot(dx, dy, color = (0.0, 0.0, 0.0), alpha=1.0, marker = ".", linestyle = "none", markersize = 1, markeredgewidth = 0.5)

            if count == 0:
                [ymin, ymax] = ax.get_ylim()
                [xmin, xmax] = ax.get_xlim()
#            ax.set_ylim(ymin, ymax)
#            ax.set_xlim(xmin, xmax)
            zoomx = xmin+(442/769)*2*xmax
            zoomy = ymin+((770-363)/770)*2*ymax
            ax.set_ylim(zoomy-50, zoomy+50)
            ax.set_xlim(zoomx-50, zoomx+50)
            plt.close()
            fig.savefig(os.path.join(Path_to_Folder,"cell"+"{}".format(count)+".png"), dpi=fig.dpi, bbox_inches='tight', pad_inches=0.0)
            
            hx.clear()
            hy.clear()
            ex.clear()
            ey.clear()
            ix.clear()
            iy.clear()
            dx.clear()
            dy.clear()

            gc.collect()                        

            index = 0
            print(count)
            count = count+1

print("Image Printing Complete")


#hc = CoordanateNested[i][n][0]
#vc = (2.0/3.0 * numpy.sin(numpy.pi/3)*(CoordanateNested[i][n][1] - CoordanateNested[i][n][2]))

