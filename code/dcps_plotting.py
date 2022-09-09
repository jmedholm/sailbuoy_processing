import matplotlib.pyplot as plt
import xarray as xr
import pandas as pd
import numpy as np
from matplotlib.dates import DateFormatter, DayLocator, HourLocator
from cmocean import cm as cmo
import matplotlib.dates as mdates

def heading_dist(ds,lr,left,right,path):

    fig,ax = plt.subplots(1,2,figsize=(20,7),sharey=True,constrained_layout=True,facecolor='w')

    for i in range(2):
        var = [left,right]
        bins = np.linspace(-lr,lr,200)
        _ = ax[i].hist(ds.where((ds['sample_heading'] > 315) | (ds['sample_heading'] < 45))[var[i]].values.flatten(), bins = bins, alpha=0.3, fc='b', edgecolor=None,label='North')
        _ = ax[i].hist(ds.where((ds['sample_heading'] > 225) & (ds['sample_heading'] < 315))[var[i]].values.flatten(), bins = bins, alpha=0.3, fc='r', edgecolor=None,label='West')
        _ = ax[i].hist(ds.where((ds['sample_heading'] > 135) & (ds['sample_heading'] < 225))[var[i]].values.flatten(), bins = bins, alpha=0.3, fc='k', edgecolor=None,label='South')
        _ = ax[i].hist(ds.where((ds['sample_heading'] > 45)  & (ds['sample_heading'] < 135))[var[i]].values.flatten(), bins = bins, alpha=0.3, fc='y', edgecolor=None,label='East')
        ax[i].legend(title='SB Direction')
        ax[i].set_title(ds[var[i]].attrs['long_name'])
        ax[i].set_xlabel(ds[var[i]].attrs['units'])
    ax[0].set_ylabel('Count')
    
    plt.savefig(path,dpi=300)


def vertical_strength_plots(ds,path):
    # Plot some random profiles of speed, std, and strength. Dim everything thats below -40dB strength
    a = 2
    b = 4
    fig, ax = plt.subplots(a,b,figsize=(30,14),constrained_layout=True,facecolor='w')
    I = np.sort(np.random.randint(len(ds.time),size=a*b))
    ax = ax.flatten()
    for i in range(len(ax)):
        axt = ax[i]
        ax1 = axt.twiny()
        (ds.isel(time=I[i])['horizontal_speed']).plot(y='depth',ylim=(83,0),alpha=1,c='C0',ax=ax1,label='Horizontal speed')
        (ds.isel(time=I[i])['sp_sd_horizontal']/100/np.sqrt(149)).plot(y='depth',ylim=(83,0),alpha=1,c='C0',ls='--',ax=ax1,label='St dev speed')
        ds.isel(time=I[i])['strength'].plot(y='depth',ylim=(83,0),alpha=1,c='C1',ax=axt)
        ax1.axhspan(ds.isel(time=I[i]).depth.where(ds.isel(time=I[i])['strength'] > -40).max()+2,83,0,1,fc='w',alpha=0.75,zorder=4)
        ax1.axhline(ds.isel(time=I[i]).depth.where(ds.isel(time=I[i])['strength'] > -40).max()+2,c='k',ls='--')
        ax1.set_xlabel('Speed (m/s)',c='C0',fontweight='bold')
        axt.set_xlabel('Strength (dB)',c='C1',fontweight='bold')
        ax1.set_title('')
        axt.set_title(str(I[i]))
        ax1.set_ylabel('')

        ax1.set_xlim(0,1)

        axt.xaxis.label.set_color('C1')        #setting up X-axis label color to yellow
        ax1.xaxis.label.set_color('C0')          #setting up Y-axis label color to blue
        axt.tick_params(axis='x', colors='C1')    #setting up X-axis tick color to red
        ax1.tick_params(axis='x', colors='C0')  #setting up Y-axis tick color to black
        ax1.spines['bottom'].set_color('C1')        # setting up Y-axis tick color to red
        ax1.spines['top'].set_color('C0')         #setting up above X-axis tick color to red
    ax1.legend(loc='lower center',framealpha=1)
    plt.savefig(path,dpi=300)

def NE_panels(ds,var,vmin,vmax,cmap,ylim,path,title,full=True):
    
    fig, ax = plt.subplots(3,1,figsize=(30,15),constrained_layout=True,sharex=True)
    
    for i in range(3):
        
        ds[var[i]].plot(y='depth',
                        ylim=ylim,
                        vmin=vmin[i], 
                        vmax=vmax[i],
                        cmap=cmap[i],
                        ax=ax[i],
                        cbar_kwargs={'label':ds[var[i]].attrs['units']})
        ax[i].set_title(ds[var[i]].attrs['description'])
        ax[i].set_ylabel('Depth (m)')
        ax[i].set_xlabel('')
    
    plt.setp(ax[0].get_xticklabels(),visible=False)
    if full == True:
        ax[1].xaxis.set_minor_locator(mdates.DayLocator([15]))
        ax[1].xaxis.set_major_locator(mdates.MonthLocator())
        ax[1].xaxis.set_major_formatter(mdates.DateFormatter("%d %b"))
        ax[1].xaxis.set_minor_formatter(mdates.DateFormatter("%d"))
    else:
        ax[1].xaxis.set_minor_locator(mdates.DayLocator())
        ax[1].xaxis.set_major_locator(mdates.DayLocator([1,5,10,15,20,25]))
        ax[1].xaxis.set_major_formatter(mdates.DateFormatter("%d %b"))
        ax[1].xaxis.set_minor_formatter(mdates.DateFormatter(""))
 
        
        
    plt.xticks(rotation=0,ha='center')
    fig.suptitle(title,fontsize='xx-large',y=1.10)

    plt.savefig(path,dpi=300)