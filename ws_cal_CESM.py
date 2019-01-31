# -*- coding: utf-8 -*-
"""
Created on Thu Jun  7 13:31:20 2018

@author: Jose Luis Rodriguez-Solis
         Oceanografia Fisica
         CICESE

WIND SPEED
Calcula los valores de wind speed basados en :

Archer, C. L., & Caldeira, K. (2008). Historical trends in the jet streams. 
Geophysical Research Letters, 35(8), 1–6. doi.org/10.1029/2008GL033614


"""

from netCDF4 import Dataset

import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

import numpy as np
#from scipy import stats

import os
import atmos

#import sys

# nivel superior
pi = 400 
# nivel infererior
ps = 100

# donde estan los archivos de Era Interim
path = '/home/luis/CESM/CESM_data/FAMIP34y2x2/cam/'
f = os.listdir(path)
f.sort()

ws = []
Ps = []
phi= []
t  = []
ls = []

nc = Dataset (path+f[0])
lat = nc.variables['lat'][:]
lon = nc.variables['lon'][:]
lev = nc.variables['lev'][:]

# buscando indices de altura de 400 y 100
zl  = np.where((lev>= 100.0) & (lev<=  400.5)); zl   = zl[0]
zs = zl[0]
zi = zl[-1]

# Malla para graficar y encontrar latitudes
x,y=np.meshgrid(lon,lat)
# buscando latitudes entre 70N y 15N
yl  = np.where((lat>= 15.0) & (lat<=  60.5)); yl   = yl[0]
ys = yl[0]
yi = yl[-1]

yeari = int(f[0][20:-15])

zl = np.zeros((len(zl),len(lat),len(lon)))
for i in range(0, len(lev[zs:zi+1])):    
    zl[i,:,:] = lev[zs+i]

num = len(f)
for i in range (0,num): # len(f) 300
    
    nc = Dataset (path+f[i])    
    u   = nc.variables['U'][:].squeeze()
    v   = nc.variables['V'][:].squeeze()
    q   = nc.variables['Q'][:].squeeze() * 1000.0
    T   = nc.variables['T'][:].squeeze() 
    nc.close()    
    
    # buscando indices de niveles superior e inferior
    us = np.nanmean( (u[:,zs:zi+1,:,:].copy())**2, axis=0)
    del u
    vs = np.nanmean( (v[:,zs:zi+1,:,:].copy())**2, axis=0)
    del v    
    ms = np.sqrt(us + vs)
    qs = np.nanmean( q[:,zs:zi+1,:,:].copy(), axis=0)
    del q 
    ts = np.nanmean( T[:,zs:zi+1,:,:].copy(), axis=0)
    del T

    # calculando rho
    rho = atmos.calculate('rho',p = zl*100.0, Tv = ts)
    
    # calculando el wind speed
    ws.append ( (np.nansum( rho * ms, axis=0)) / (np.nansum(rho, axis=0) ) )
    # calculando Presion
    Ps.append ( (np.nansum( rho * ms * zl, axis=0 )) / (np.nansum(rho * ms, axis=0) ) )
    # calculando latitud
    ls_top = np.nansum ( rho * ms, axis=0) * y
    ls_top = np.nansum ( ls_top[ys:yi,:], axis = 0 )
    
    ls_inf = np.nansum ( rho * ms, axis=0)
    ls_inf = np.nansum ( ls_inf[ys:yi,:], axis = 0 )
    
    ls.append( ls_top/ls_inf )
    
    t.append(f[i][7:-3])
    print f[i][7:-3]

#yearf = int(f[i][7:-7])
yearf = int(f[0][20:-15])

print ('finaliza ciclo de calculos')
print ('')
print ('comienza ciclo para medias')
#%% inicia promedios 

ws = np.asarray(ws)
wsm = np.nanmean(ws,axis = 0)

wn = np.nanmean (ws[:,ys:yi,:], axis = 1)
wn = np.nanmean (wn, axis = 1)

Ps = np.asarray(Ps)
psm = np.nanmean(Ps,axis = 0)

pn = np.nanmean (Ps[:,ys:yi,:], axis = 1)
pn = np.nanmean (pn, axis = 1)

ls = np.asarray(ls)
ls = np.nanmean(ls,axis=1)
#%% VIENTO
# promedios para verano
print ('viento')
wn_ver = []
wnsver = []
k = 0
for j in range(5,len(pn),12):
    print t[j:j+3]
    wn_ver.append(np.nanmean (wn[j:j+3]))
    wnsver.append(np.nanmean (ws[j:j+3,:,:],axis=0 ) )
wnsver = np.asarray( wnsver )
# promedios para invierno
wn_inv = []
wnsinv = []
k = 0
for j in range(11,len(wn)-1,12):
    print t[j:j+3]
    wn_inv.append(np.nanmean (wn[j:j+3]))    
    wnsinv.append(np.nanmean (ws[j:j+3,:,:],axis=0 ) )
wnsinv = np.asarray( wnsinv )
# promedio anual    
wni = []
k = 0
for i in range(12,num+1,12):
    print i
    wni.append( np.nanmean(ws[k:i,ys:yi+1,:] ) )
    k = k + 12
wni = np.asarray(wni)
    
#%% PRESION
# promedios para verano
print ('presion')
pn_ver = []
psver  = []
k = 0
for j in range(5,len(pn),12):
    print t[j:j+3]
    pn_ver.append(np.nanmean (pn[j:j+3]))
    psver.append ( np.nanmean(Ps[j:j+3,:,:],axis = 0) )
psver = np.asarray(psver)
# promedios para invierno
pn_inv = []
psinv  = []
k = 0
for j in range(11,len(pn)-1,12):
    print t[j:j+3]
    pn_inv.append(np.nanmean (pn[j:j+3]))    
    psinv.append ( np.nanmean(Ps[j:j+3,:,:],axis = 0) )
psinv = np.asarray(psinv)
# promedio anual    
pni = []
k = 0
for i in range(12,num+1,12):
    print i
    pni.append( np.nanmean(Ps[k:i,ys:yi+1,:] ) )
    k = k + 12
pni = np.asarray(pni)    
    
#%% LATITUD
print ('Latitud')
# promedios para latitud anual
ln = []
k = 0
for j in range(0,len(pn),12):
    print t[j:j+12]
    ln.append(np.nanmean (ls[j:j+12]))    
ln = np.asarray(ln)
# promedios para verano
ln_ver = []
k = 0
for j in range(5,len(pn),12):
    print t[j:j+3]
    ln_ver.append(np.nanmean (ls[j:j+3]))
ln_ver = np.asarray( ln_ver )
# promedios para invierno
ln_inv = []
k = 0
for j in range(11,len(wn)-1,12):
    print t[j:j+3]
    ln_inv.append(np.nanmean (ls[j:j+3]))    
ln_inv = np.asarray( ln_inv )
#%%
ya = -90           # latitud inicial
yf =  90            # latitud final
xa = 0           # longitud inicial
xf = 360           # longitud final

#%% promedio anual viento
print ('comienza a graficar')
fig, ax = plt.subplots(num=1, figsize=(8,4))
m = Basemap(projection='cyl',llcrnrlat=ya,urcrnrlat=yf,\
            llcrnrlon=xa,urcrnrlon=xf,resolution='l')
m.drawcoastlines(linewidth=0.5)
m.drawcountries(linewidth=0.5)
parallels = np.arange(ya,yf+0.1,45.)
m.drawparallels(parallels,labels=[1,0,0,0],linewidth= 0.2,fontsize=8)
meridians = np.arange(xa,xf,45.)
m.drawmeridians(meridians,labels=[0,0,0,1],linewidth= 0.2,fontsize=8)
clev = np.arange(10,46,5)
cs = plt.contourf(x,y,wsm,shading='interp',levels=clev,cmap=plt.cm.Greens,latlon=True,extend='both')
cbar = m.colorbar(cs,location='bottom',ticks=clev,size="05%",pad="15%",extend='both',drawedges=True)
cbar.set_label(u'm/s',x=1.05,fontsize=8)
plt.title ('Promedio anual de velocidad del viento (Jet Stream)\nDatos:CESM-AMIP 1981-2005',fontsize=8,loc = 'left')
plt.savefig('CAMcesmAnnualMeanWS.png',bbox_inches='tight',dpi=150)
plt.close()

fig, ax = plt.subplots(num=2, figsize=(8,4))
m = Basemap(projection='cyl',llcrnrlat=ya,urcrnrlat=yf,\
            llcrnrlon=xa,urcrnrlon=xf,resolution='l')
m.drawcoastlines(linewidth=0.5)
m.drawcountries(linewidth=0.5)
parallels = np.arange(ya,yf+0.1,45.)
m.drawparallels(parallels,labels=[1,0,0,0],linewidth= 0.2,fontsize=8)
meridians = np.arange(xa,xf,45.)
m.drawmeridians(meridians,labels=[0,0,0,1],linewidth= 0.2,fontsize=8)
clev = np.arange(10,46,5)
cs = plt.contourf(x,y,np.nanmean(wnsinv,axis=0),shading='interp',levels=clev,cmap=plt.cm.Greens,latlon=True,extend='both')
cbar = m.colorbar(cs,location='bottom',ticks=clev,size="05%",pad="15%",extend='both',drawedges=True)
cbar.set_label(u'm/s',x=1.05,fontsize=8)
plt.title ('Promedio de velocidad del viento en Invierno (Jet Stream)\nDatos:CESM-AMIP 1981-2005',fontsize=8,loc = 'left')
plt.savefig('CAMWinterMeanWS.png',bbox_inches='tight',dpi=150)
plt.close()

fig, ax = plt.subplots(num=3, figsize=(8,4))
m = Basemap(projection='cyl',llcrnrlat=ya,urcrnrlat=yf,\
            llcrnrlon=xa,urcrnrlon=xf,resolution='l')
m.drawcoastlines(linewidth=0.5)
m.drawcountries(linewidth=0.5)
parallels = np.arange(ya,yf+0.1,45.)
m.drawparallels(parallels,labels=[1,0,0,0],linewidth= 0.2,fontsize=8)
meridians = np.arange(xa,xf,45.)
m.drawmeridians(meridians,labels=[0,0,0,1],linewidth= 0.2,fontsize=8)
clev = np.arange(10,46,5)
cs = plt.contourf(x,y,np.nanmean(wnsver,axis=0),shading='interp',levels=clev,cmap=plt.cm.Greens,latlon=True,extend='both')
cbar = m.colorbar(cs,location='bottom',ticks=clev,size="05%",pad="15%",extend='both',drawedges=True)
cbar.set_label(u'm/s',x=1.05,fontsize=8)
plt.title ('Promedio de velocidad del viento en Verano (Jet Stream)\nDatos:CESM-AMIP 1981-2005',fontsize=8,loc = 'left')
plt.savefig('CAMSummerMeanWS.png',bbox_inches='tight',dpi=150)
plt.close()
#%% promedio anual presion
fig, ax = plt.subplots(num=4, figsize=(8,4))
m = Basemap(projection='cyl',llcrnrlat=ya,urcrnrlat=yf,\
            llcrnrlon=xa,urcrnrlon=xf,resolution='l')
m.drawcoastlines(linewidth=0.5)
m.drawcountries(linewidth=0.5)
parallels = np.arange(ya,yf+0.1,45.)
m.drawparallels(parallels,labels=[1,0,0,0],linewidth= 0.2,fontsize=8)
meridians = np.arange(xa,xf,45.)
m.drawmeridians(meridians,labels=[0,0,0,1],linewidth= 0.2,fontsize=8)
clev = np.arange(230,301,5)
cs = plt.contourf(x,y,psm,shading='interp',levels=clev,cmap=plt.cm.bone_r,latlon=True)
cbar = m.colorbar(cs,location='bottom',ticks=np.arange(230,301,10),size="05%",pad="15%",extend='both',drawedges=True)
cbar.set_label(u'hPa', x=1.05,fontsize=8)
plt.title (u'Promedio anual de Presión (Jet Stream)\nDatos:CESM-AMIP 1981-2005', fontsize=8,loc = 'left')
plt.savefig('CAMAnnualMeanP.png',bbox_inches='tight',dpi=150)
#plt.close()

fig, ax = plt.subplots(num=5, figsize=(8,4))
m = Basemap(projection='cyl',llcrnrlat=ya,urcrnrlat=yf,\
            llcrnrlon=xa,urcrnrlon=xf,resolution='l')
m.drawcoastlines(linewidth=0.5)
m.drawcountries(linewidth=0.5)
parallels = np.arange(ya,yf+0.1,45.)
m.drawparallels(parallels,labels=[1,0,0,0],linewidth= 0.2,fontsize=8)
meridians = np.arange(xa,xf,45.)
m.drawmeridians(meridians,labels=[0,0,0,1],linewidth= 0.2,fontsize=8)
clev = np.arange(230,301,5)
cs = plt.contourf(x,y,np.nanmean(psinv,axis=0),shading='interp',levels=clev,cmap=plt.cm.bone_r,latlon=True)
cbar = m.colorbar(cs,location='bottom',ticks=np.arange(230,301,10),size="05%",pad="15%",extend='both',drawedges=True)
cbar.set_label(u'hPa', x=1.05,fontsize=8)
plt.title (u'Promedio de Presión de Invierno (Jet Stream)\nDatos:CESM-AMIP 1981-2005', fontsize=8,loc = 'left')
plt.savefig('CAMWinterMeanP.png',bbox_inches='tight',dpi=150)
#plt.close()

fig, ax = plt.subplots(num=6, figsize=(8,4))
m = Basemap(projection='cyl',llcrnrlat=ya,urcrnrlat=yf,\
            llcrnrlon=xa,urcrnrlon=xf,resolution='l')
m.drawcoastlines(linewidth=0.5)
m.drawcountries(linewidth=0.5)
parallels = np.arange(ya,yf+0.1,45.)
m.drawparallels(parallels,labels=[1,0,0,0],linewidth= 0.2,fontsize=8)
meridians = np.arange(xa,xf,45.)
m.drawmeridians(meridians,labels=[0,0,0,1],linewidth= 0.2,fontsize=8)
clev = np.arange(230,301,5)
cs = plt.contourf(x,y,np.nanmean(psver,axis=0),shading='interp',levels=clev,cmap=plt.cm.bone_r,latlon=True)
cbar = m.colorbar(cs,location='bottom',ticks=np.arange(230,301,10),size="05%",pad="15%",extend='both',drawedges=True)
cbar.set_label(u'hPa', x=1.05,fontsize=8)
plt.title (u'Promedio de Presión de Verano (Jet Stream)\nDatos:CESM-AMIP 1981-2005', fontsize=8,loc = 'left')
plt.savefig('CAMSummerMeanP.png',bbox_inches='tight',dpi=150)
#plt.close()

#%% serie de tendencia de Presion 
ti  = np.arange(1981,yearf+1,1/12.0)
tip = np.arange(1981,yearf+1,1)
fig, ax = plt.subplots(num=7, figsize=(6,8))

yt = np.arange(0,len(tip),1)

plt.subplot(211)

za = np.polyfit(yt,pni,1)
zv = np.polyfit(yt,pn_ver,1)
zi = np.polyfit(yt[0:-1],pn_inv,1)

#plt.plot(ti,pn,'k',linewidth=0.5)
plt.plot(tip,pni,'k',linewidth=0.5)
plt.plot(tip,(za[0]*yt+za[1]),':k',linewidth=1.0)
plt.plot(tip,pn_ver,'r',linewidth=0.5)
plt.plot(tip,(zv[0]*yt+zv[1]),':r',linewidth=1.0)
plt.plot(tip[0:-1],pn_inv,'b',linewidth=0.5)
plt.plot(tip[0:-1],(zi[0]*yt[0:-1]+zi[1]),':b',linewidth=1.0)
#plt.plot([yeari,yearf+1],[np.nanmean(pn),np.nanmean(pn)], ':k',linewidth=0.5)
plt.ylim(255,242)
plt.xlim(yeari,yearf+1)
plt.ylabel(u'Presión [hPa]')
plt.title (u'Comportamiento de presión y viento de Jet Stream\nDatos:CESM-AMIP 1981-2017', fontsize=8,loc = 'left')
plt.text(2010,241.0,'-Anual',fontsize=8)
plt.text(2010,241.5,'-Verano',fontsize=8,color='r')
plt.text(2010,242,'-Invierno',fontsize=8,color='b')

plt.subplot(212)

za = np.polyfit(yt,wni,1)
zi = np.polyfit(yt[0:-1],wn_inv,1)
zv = np.polyfit(yt,wn_ver,1)

plt.plot(tip,wni,'k',linewidth=0.5)
plt.plot(tip,(za[0]*yt+za[1]),':k',linewidth=1.0)
plt.plot(tip,wn_ver,'r',linewidth=0.5)
plt.plot(tip,(zv[0]*yt+zv[1]),':r',linewidth=1.0)
plt.plot(tip[0:-1],wn_inv,'b',linewidth=0.5)
plt.plot(tip[0:-1],(zi[0]*yt[0:-1]+zi[1]),':b',linewidth=1.0)
#plt.plot([yeari,yearf+1],[np.nanmean(wn),np.nanmean(wn)], ':k',linewidth=0.5)
recta = 'y='+str(np.round(za[0],4))+'x+'+str(np.round(za[1],4))
plt.text(2008,za[1]+1.5,recta,fontsize=6,fontweight='bold')
recta = 'y='+str(np.round(zv[0],4))+'x+'+str(np.round(zv[1],4))
plt.text(2008,zv[1]+1.5,recta,fontsize=6,fontweight='bold')
recta = 'y='+str(np.round(zi[0],4))+'x+'+str(np.round(zi[1],4))
plt.text(2008,zi[1]+1.5,recta,fontsize=6,fontweight='bold')
plt.ylim(10,30)
plt.xlim(yeari,yearf+1)
plt.ylabel('Viento [m/s]')
plt.xlabel('Tiempo')

plt.savefig('CAMMeanSeries.png',bbox_inches='tight',dpi=150)
#%%
# serie de tendencia de latitud
fig, ax = plt.subplots(num=8, figsize=(6,8))

plt.subplot(311)

za = np.polyfit(yt,ln-np.nanmean(ln),1)
zi = np.polyfit(yt[0:-1],ln_inv-np.nanmean(ln_inv),1)
zv = np.polyfit(yt,ln_ver-np.nanmean(ln_ver),1)

plt.plot(tip,ln-np.nanmean(ln),'k',linewidth=0.5)
plt.plot(tip,(za[0]*yt+za[1]),':k',linewidth=1.0)
plt.ylim(-1.9,1.9)
plt.ylabel(u'Lat Anomalía [deg]')
plt.text(1982,1.6,u'Anual',fontsize=10)
recta = 'y='+str(np.round(za[0],4))+'x'+str(np.round(za[1],4))
plt.text(2008,1.6,recta,fontsize=8,fontweight='bold')
plt.title (u'Anomalía de la posición del Jet Stream\nDatos:CESM-AMIP 1981-2017', fontsize=8,loc = 'left')
plt.xlim(yeari,yearf+1)

plt.subplot(312)
plt.plot(tip,ln_ver-np.nanmean(ln_ver),'r',linewidth=0.5)
plt.plot(tip,(zv[0]*yt+zv[1]),':r',linewidth=1.0)
plt.ylim(-1.9,1.9)
plt.ylabel(u'Lat Anomalía [deg]')
plt.text(1982,1.6,u'Verano',fontsize=10)
recta = 'y='+str(np.round(zv[0],4))+'x'+str(np.round(zv[1],4))
plt.text(2008,1.6,recta,fontsize=8,fontweight='bold')
plt.xlim(yeari,yearf+1)

plt.subplot(313)
plt.plot(tip[0:-1],ln_inv-np.nanmean(ln_inv),'b',linewidth=0.5)
plt.plot(tip[0:-1],(zi[0]*yt[0:-1]+zi[1]),':b',linewidth=1.0)
#plt.plot([yeari,yearf+1],[0,0], ':k',linewidth=0.5)
plt.ylim(-1.9,1.9)
plt.xlim(yeari,yearf+1)
plt.ylabel(u'Lat Anomalía [deg]')
plt.xlabel('tiempo')
plt.text(1982,1.6,u'Invierno',fontsize=10)
recta = 'y='+str(np.round(zi[0],4))+'x'+str(np.round(zi[1],4))
plt.text(2008,1.6,recta,fontsize=8,fontweight='bold')
plt.savefig('CAMLatMeanSeries.png',bbox_inches='tight',dpi=150)
