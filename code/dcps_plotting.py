import matplotlib.pyplot as plt
import xarray as xr
import pandas as pd
import numpy as np
from matplotlib.dates import DateFormatter, DayLocator, HourLocator
from cmocean import cm as cmo
import matplotlib.dates as mdates
from matplotlib.gridspec import GridSpec
import cartopy.crs as ccrs
import matplotlib
import gsw

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
    fig, ax = plt.subplots(a,b,figsize=(30,10),constrained_layout=True,facecolor='w',sharex='col',sharey='row')
    I = np.sort(np.random.randint(len(ds.time),size=a*b))
    ax = ax.flatten()
    for i in range(len(ax)):
        axt = ax[i]
        ax1 = axt.twiny()
        (ds.isel(time=I[i])['horizontal_speed']).plot(y='depth',ylim=(83,0),alpha=1,c='C0',ax=ax1,label='Horizontal speed')
        (ds.isel(time=I[i])['sp_sd_horizontal']/100/np.sqrt(149)).plot(y='depth',ylim=(83,0),alpha=1,c='C0',ls='--',ax=ax1,label='St dev speed')
        ds.isel(time=I[i])['strength'].plot(y='depth',ylim=(83,0),alpha=1,c='C1',ax=axt)
        ax1.axhspan(ds.isel(time=I[i]).depth.where(ds.isel(time=I[i])['strength'] > -44).max()+2,83,0,1,fc='w',alpha=0.75,zorder=4)
        ax1.axhline(ds.isel(time=I[i]).depth.where(ds.isel(time=I[i])['strength'] > -44).max()+2,c='k',ls='--')
        ax1.set_xlabel('Speed (m/s)',c='C0',fontweight='bold',labelpad=15)
        axt.set_xlabel('Strength (dB)',c='C1',fontweight='bold')
        ax1.set_title('')
        axt.set_title('')
        axt.set_title(f"{I[i]}  ",loc='right',y=0.85)
        ax1.set_ylabel('')
        ax1.set_xlim(0,1)
        axt.xaxis.label.set_color('C1')        #setting up X-axis label color to yellow
        ax1.xaxis.label.set_color('C0')          #setting up Y-axis label color to blue
        axt.tick_params(axis='x', colors='C1')    #setting up X-axis tick color to red
        ax1.tick_params(axis='x', colors='C0')  #setting up Y-axis tick color to black
        ax1.spines['bottom'].set_color('C1')        # setting up Y-axis tick color to red
        ax1.spines['top'].set_color('C0')         #setting up above X-axis tick color to red
        if i in [0,1,2,3]:
            axt.set_xlabel('')
            
        if i in [4,5,6,7]:
            ax1.set_xlabel('')
            plt.setp(ax1.get_xticklabels(), visible=False)
        if i in [0,4]:
            axt.set_ylabel('Depth (m)')
        else:
            axt.set_ylabel('')
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

def rot_ticks(axs,rot,ha):
    for xlabels in axs.get_xticklabels():
                xlabels.set_rotation(rot)
                xlabels.set_ha(ha)

def plot_bathymetry(ax,bathy,color=True):
    
    ct = bathy.plot.contour(levels=25,
                       linestyles='-',
                       linewidths=1,
                       colors='k',
                       ax=ax,
                       add_colorbar=False,
                       transform=ccrs.PlateCarree(),
                       zorder=2,
                       vmin=1000,
                       vmax=6000)
    if color == True:
        ct = bathy.plot.contourf(levels=25,
                                 ax=ax,
                                 cmap='Blues',
                                 add_colorbar=False,
                                 transform=ccrs.PlateCarree(),
                                 zorder=1,
                                 vmin=1000,
                                 vmax=6000)

    
    return ct

def plot_subset_map(ds,path,bathy,xlims,ylims):
    fig = plt.figure(figsize=(12,15),facecolor='w')

    spec = GridSpec(ncols=2,
                    nrows=20,
                    figure=fig,
                    wspace=0.05,
                    hspace=0.5)

    ax = fig.add_subplot(spec[:18, :],
                         projection=ccrs.SouthPolarStereo())

    cax1 = fig.add_subplot(spec[19, 0])
    cax2 = fig.add_subplot(spec[19, 1])

    ct = plot_bathymetry(ax,bathy)

    sct = ax.scatter(ds.longitude,
                     ds.latitude,
                     c=mdates.date2num(ds.time.values),
                     cmap='spring',
                     transform=ccrs.PlateCarree(),
                     zorder=3)

    cb1 = plt.colorbar(sct,
                       cax=cax1,
                       orientation='horizontal',
                       pad=0.01,
                       aspect=50,
                       fraction=0.025)

    cb1.ax.xaxis.set_major_locator(mdates.DayLocator([1]))
    cb1.ax.xaxis.set_minor_locator(mdates.DayLocator([10,20]))
    cb1.ax.xaxis.set_major_formatter(mdates.DateFormatter("%d %b"))

    cb2 = plt.colorbar(ct,
                       cax=cax2,
                       orientation='horizontal',
                       pad=0.01,
                       aspect=50,
                       fraction=0.025,
                       label='Depth (m)')

    cb2.set_ticks([2500,5000])


    ax.set_extent([xlims[0], xlims[1], ylims[0], ylims[1]], ccrs.PlateCarree())

    gl  = ax.gridlines(draw_labels=True,
                       dms=False, 
                       x_inline=False, 
                       y_inline=False,
                       color='k',
                       linestyle='--',
                       transform=ccrs.PlateCarree(),
                       zorder=2)

    gl.xlines          = True
    gl.ylines          = True
    gl.top_labels      = False
    gl.bottom_labels   = True
    gl.left_labels     = True
    gl.right_labels    = False

    cax1.tick_params(axis="x", which="major", length=6,width=1)
    cax2.tick_params(axis="x", which="major", length=6,width=1)
    cb2.ax.invert_xaxis()
    cb2.ax.xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(500))
    title = f"SB Kringla, {ds['time'][0].values.astype('datetime64[D]')} to {ds['time'][-1].values.astype('datetime64[D]')}"
    ax.set_title('SB Kringla',pad=10)

    plt.savefig(f"{path}Subset_map.png")
    
def panel_map_plot_time(ds,dates,i,c,l,bathy,xlims,ylims,path):
    
    fig = plt.figure(figsize=(30,10),facecolor='w')
    gs = fig.add_gridspec(nrows=6, ncols=40, hspace=0.3, wspace=0)
    
    axt = fig.add_subplot(gs[0:2,:28])
    axs = axt.twinx()
    ax1 = fig.add_subplot(gs[2:4, :28],sharex=axt)
    ax2 = fig.add_subplot(gs[4:, :28],sharex=axt)
    ax3 = fig.add_subplot(gs[1:5, 33:])
    cax1 = fig.add_subplot(gs[2:4, 29])
    cax2 = fig.add_subplot(gs[4:, 29])
    
    ds.sel(time=slice(dates[i],dates[i+1]))['sst'].plot(x='time',ax=axt,c='C1')
    ds.sel(time=slice(dates[i],dates[i+1]))['sss'].plot(x='time',ax=axs,c='C0')
    
    ds['vel_north'].sel(time=slice(dates[i],dates[i+1])).plot(x='time',
                                                                     y='depth',
                                                                     ylim=(59,0),
                                                                     vmin=-0.75,
                                                                     vmax=0.75,
                                                                     ax=ax1,
                                                                     cmap='cmo.balance',
                                                                     cbar_kwargs={'pad':0.01,
                                                                                  'label':'North speed (m s$^{-1}$)',
                                                                                  'cax':cax1})
    ds['vel_east'].sel(time=slice(dates[i],dates[i+1])).plot(x='time',
                                                                 y='depth',
                                                                 ylim=(59,0),
                                                                 vmin=-0.75,
                                                                 vmax=0.75,
                                                                 ax=ax2,
                                                                 cmap='cmo.balance',
                                                                 cbar_kwargs={'pad':0.01,
                                                                              'label':'East speed (m s$^{-1}$)',
                                                                              'cax':cax2})
    plt.setp(ax1.get_xticklabels(),visible=False)
    plt.setp(axt.get_xticklabels(),visible=False)
    ax1.set_xlabel('')
    ax1.set_ylabel('Depth (m)')
    ax2.set_ylabel('Depth (m)')
    #ax1.set_xlim(-150,150)
    #ax1.set_xlim(-150,150)    
    ax3.scatter(ds.longitude.sel(time=slice(dates[i],dates[i+1])),
                ds.latitude.sel(time=slice(dates[i],dates[i+1])),
                c=ds.latitude.sel(time=slice(dates[i],dates[i+1])).time,
                zorder=4,
                label=l[i])
    
    bathy.plot.contour(levels=[3500],colors='k',ax=ax3,zorder=1)
    bathy.plot.contour(colors='gray',linewidths=1,ax=ax3,zorder=1)
    ax3.set_xlim(xlims+np.array([0,1]))
    ax3.set_ylim(ylims)
    ax3.set_xlabel('Longitude (°E)')
    ax3.set_ylabel('Latitude (°N)',labelpad=10)
    ax3.yaxis.tick_right()
    ax3.yaxis.set_label_position("right")
    
    axt.set_title(l[i])
    
    dur = ds.sel(time=slice(dates[i],dates[i+1])).time.diff('time').sum().values.astype('timedelta64[h]')
    mid = ds.sel(time=slice(dates[i],dates[i+1])).time[0] + dur/2
    xdates = [(mid-np.timedelta64(85,'h')).values,(mid+np.timedelta64(85,'h')).values]
    
    rot_ticks(ax2,0,'center')
    ax3.set_title(dur)
    ax2.set_xlim(xdates)
    ax1.xaxis.set_major_locator(mdates.HourLocator([0]))
    ax1.xaxis.set_minor_locator(mdates.HourLocator([0,3,6,9,12,15,18,21]))
    ax1.xaxis.set_major_formatter(mdates.DateFormatter("%d %b\n%H:%M"))
    ax1.xaxis.set_minor_formatter(mdates.DateFormatter(""))
    ax1.set_xlabel('')
    axt.set_xlabel('')
    ax2.set_xlabel('')
    axt.set_ylim(-0.5,1.25)
    axs.set_ylim(33.5,34.2)
    axs.set_ylabel('Salinity (g kg$^{-1}$)', c='C0')
    axt.set_ylabel('Temperature (°C)', c='C1')
    axs.tick_params(axis='y', color='C0', labelcolor='C0')
    axt.tick_params(axis='y', color='C1', labelcolor='C1')
    axs.spines['left'].set_edgecolor('C1')
    axs.spines['right'].set_edgecolor('C0')
    axs.spines['left'].set_linewidth(2)
    axs.spines['right'].set_linewidth(2)
    
    plt.savefig(f"{path}_{c:02d}.png")
    
def distance_plot(ds_across,i,bathy,xlims,ylims,path,name,title):
        
    fig = plt.figure(figsize=(30,10),facecolor='w')
    gs = fig.add_gridspec(nrows=6, ncols=40, hspace=0.3, wspace=0)
    
    axt = fig.add_subplot(gs[0:2,:28])
    axs = axt.twinx()
    ax1 = fig.add_subplot(gs[2:4, :28],sharex=axt)
    ax2 = fig.add_subplot(gs[4:, :28],sharex=axt)
    ax3 = fig.add_subplot(gs[1:5, 33:])
    cax1 = fig.add_subplot(gs[2:4, 29])
    cax2 = fig.add_subplot(gs[4:, 29])
    
    ds_across[i]['sst'].plot(x='distance',ax=axt,c='C1')
    ds_across[i]['sss'].plot(x='distance',ax=axs,c='C0')
    ds_across[i]['vel_north'].plot(y='depth',ylim=(59,0),vmin=-0.75,vmax=0.75,ax=ax1,cmap='cmo.balance',cbar_kwargs={'pad':0.01,'label':'North speed (m s$^{-1}$)','cax':cax1})
    ds_across[i]['vel_east'].plot(y='depth',ylim=(59,0),vmin=-0.75,vmax=0.75,ax=ax2,cmap='cmo.balance',cbar_kwargs={'pad':0.01,'label':'East speed (m s$^{-1}$)','cax':cax2})
    plt.setp(ax1.get_xticklabels(),visible=False)
    plt.setp(axt.get_xticklabels(),visible=False)
    axt.set_xlabel('')
    ax1.set_xlabel('')
    ax1.set_ylabel('Depth (m)')
    ax2.set_ylabel('Depth (m)')
    axt.set_xlim(150,-150)
    ax1.set_xlim(150,-150)
    ax2.set_xlim(150,-150)
    axt.set_ylim(-0.25,1.25)
    axs.set_ylim(33.5,34.2)
    axs.set_ylabel('Salinity (g kg$^{-1}$)', c='C0')
    axt.set_ylabel('Temperature (°C)', c='C1')
    axs.tick_params(axis='y', color='C0', labelcolor='C0')
    axt.tick_params(axis='y', color='C1', labelcolor='C1')
    axs.spines['left'].set_edgecolor('C1')
    axs.spines['right'].set_edgecolor('C0')
    axs.spines['left'].set_linewidth(2)
    axs.spines['right'].set_linewidth(2)

    xr.plot.scatter(ds_across[i],'longitude','latitude',c=ds_across[i].time,zorder=2,ax=ax3)
    bathy.plot.contour(levels=[3500],colors='k',ax=ax3,zorder=1)
    bathy.plot.contour(colors='gray',linewidths=1,ax=ax3,zorder=1)

    ax3.set_xlim(xlims+np.array([0,1]))
    ax3.set_ylim(ylims)
    ax3.set_xlabel('Longitude (°E)')
    ax3.set_ylabel('Latitude (°N)',labelpad=10)
    ax3.yaxis.tick_right()
    ax3.yaxis.set_label_position("right")
    #ax3.tick_params(axis='both', which='major', pad=0)

    dur = (ds_across[i].time.values[-1].astype('datetime64[h]') - ds_across[i].time.values[0].astype('datetime64[h]'))
    axt.set_title(f"Across-slope {i+1}, {title}")
    ax3.set_title(dur)
    plt.savefig(f"{path}/across-slope/gridded_distance/across_{(i+1):02d}_{name}.png")
    
def full_plots(ds,path,name,vlim=0.1,hlim=0.0005):
    
    for j in range(len(ds)):
        
        fig, ax = plt.subplots(8,1,figsize=(30,21),sharex=True,constrained_layout=True)

        var = ['wind_speed','air_temp','sst','vel_north','vel_east','vert_shear','hor_shear','hor_shear2']
        ylims = [[0,20],[-5,5],[-1,1],[50,0],[50,0],[50,0],[50,0],[50,0]]
        clims = [[],[],[],[-0.75,0.75],[-0.75,0.75],[0,vlim],[0,hlim],[0,hlim]]
        cmaps = [[],[],[],cmo.balance,cmo.balance,cmo.tempo,cmo.tempo,cmo.tempo]
        cbar_fmt = ['','','','%.2f','%.2f','%.3f','%.5f','%.5f']
        for i, axs in enumerate(ax):
            if i < 3:
                props = {'ylim':ylims[i], 'lw':3, 'color':'C1'}
                ds[j][var[i]].sortby('distance').plot(x='distance',xlim=(150,-150),ax=axs,**props)
                axs.set_ylabel(f"({ds[j][var[i]].attrs['units']})",c=['k','k','C1'][i])

            else:
                props = {'ylim':ylims[i], 'vmin':clims[i][0], 'vmax':clims[i][1], 'cmap':cmaps[i]}
                cbar_kwargs={'label':f"({ds[j][var[i]].attrs['units']})",'pad':-0.035,'format':cbar_fmt[i],'aspect':10}
                ds[j][var[i]].sortby('distance').plot(y='depth',x='distance',xlim=(150,-150),ax=axs,cbar_kwargs=cbar_kwargs,**props)
                axs.set_ylabel(f"Depth (m)")

            axs.set_xlabel('')

            if i == 2:
                axs.set_title(f"({['a','b','c','d','e','f','g','h'][i]}) {ds[j][var[i]].attrs['long_name']} and salinity",loc='left',y=.05,x=0.01)
            else:
                axs.set_title(f"({['a','b','c','d','e','f','g','h'][i]}) {ds[j][var[i]].attrs['long_name']}",loc='left',y=.05,x=0.01)
            if i == 0:
                axs.set_title(f"({['a','b','c','d','e','f','g','h'][i]}) {ds[j][var[i]].attrs['long_name']} and gusts",loc='left',y=.05,x=0.01)

        ax_s = ax[2].twinx()
        ds[j]['sss'].sortby('distance').plot(x='distance',xlim=(150,-150),ax=ax_s,c='C0',lw=3)
        ax_s.set_ylabel(f"({ds[j]['sss'].attrs['units']})",c='C0')

        for i in range(2):
            [ax[2],ax_s][i].tick_params(axis='y', color=['C1','C0'][i], labelcolor=['C1','C0'][i])
            ax_s.spines[['left','right'][i]].set_edgecolor(['C1','C0'][i])
            ax_s.spines[['left','right'][i]].set_linewidth(2)

        a_b = gsw.alpha_on_beta(ds[j]['sss'].mean('time'),ds[j]['sst'].mean('time'),0).values
        ax_s.set_ylim(33.3,34.2)
        ax[2].set_ylim(-2,-2+(1/a_b)*0.9)
        ax[0].fill_between(ds[j]['distance'],ds[j]['wind_speed'],ds[j]['wind_gust'],fc='C1',ec=None,alpha=0.25)

        axs.set_xlabel('Distance (km)')
        fig.suptitle(f"Across-slope {[1,2,3,4,5,6,7,8,9,10,11,12,15,16][j]}",fontsize='xx-large')
        plt.savefig(f"{path}/across-slope/panels_across_{[1,2,3,4,5,6,7,8,9,10,11,12,15,16][j]:02d}_{name}.png")