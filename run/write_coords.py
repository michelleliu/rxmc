import numpy as np
import random

random.seed()

Npart=700
Nside=0
while (Nside*Nside*Nside<Npart):
    Nside+=1
print "Npart=%d, Nside=%d" % (Npart,Nside)
Ntype=2
Nbox=2
Boxsize=[10, 10]

Npboxtype=np.zeros([Nbox,Ntype],dtype=int)
Npboxtype[0,0]=130
Npboxtype[1,0]=70
Npboxtype[0,1]=370
Npboxtype[1,1]=130
Box=np.zeros(Nbox,dtype=float)

Rx=np.zeros(Npart,dtype=float)
Ry=np.zeros(Npart,dtype=float)
Rz=np.zeros(Npart,dtype=float)

Types=np.zeros(Npart,dtype=int)
Ibox=np.zeros(Npart,dtype=int)

for i in range(Nbox):
    Box[i]=Boxsize[i]

xarray=np.zeros(Nside*Nside*Nside,dtype=int)
yarray=np.zeros(Nside*Nside*Nside,dtype=int)
zarray=np.zeros(Nside*Nside*Nside,dtype=int)
latticecounter=-1
for i in range(Nside):
    for j in range(Nside):
        for k in range(Nside):
            latticecounter+=1
            xarray[latticecounter]=i
            yarray[latticecounter]=j
            zarray[latticecounter]=k

Coordold=open("Coordold","w+")
Coordold.write("{0} {1}\n".format(Box[0],Box[1]))
Coordold.write("{0}\n".format(Npart))
Coordold.write("{0}\n".format(Ntype))

icounter=-1
for typeidx in range(Ntype):
    for boxidx in range(Nbox):
        for ipart in range(Npboxtype[boxidx,typeidx]):
            icounter+=1
            Rx[icounter]=xarray[icounter]*float(Boxsize[boxidx])/float(Nside)
            Ry[icounter]=yarray[icounter]*float(Boxsize[boxidx])/float(Nside)
            Rz[icounter]=zarray[icounter]*float(Boxsize[boxidx])/float(Nside)
            Types[icounter]=typeidx+1
            Ibox[icounter]=boxidx+1
            #print "{0} {1} {2} {3} {4}\n".format(Rx[icounter],Ry[icounter],
            #    Rz[icounter],Types[icounter],Ibox[icounter])
            Coordold.write("{0} {1} {2} {3} {4}\n".format(Rx[icounter],Ry[icounter],
                Rz[icounter],Types[icounter],Ibox[icounter]))
