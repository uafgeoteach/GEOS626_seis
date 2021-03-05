import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pylab as plt
def plot_event_station(elat,elon,w=[],slat=[],slon=[],stas=[]):
    '''
    Input:
        elat:   event latitude
        elon:   event longitude
        w:      stream containing station information
        slat:   if not given w, array of station latitudes to plot
        slon:   if not given w, array of station longitudes to plot
        stas:   if not given w, array of station names to plot
    '''
    
    # Plot global map of selected stations
    lat_epic = elat    # Latitude
    lon_epic = elon   # Longitude
    
    
    fig = plt.figure(figsize=[15, 8])
    ax = fig.add_subplot(1, 1, 1,
                         projection=ccrs.AzimuthalEquidistant(central_longitude=lon_epic,
                                                              central_latitude=lat_epic))
    
    ax.set_global()
    ax.coastlines()
    ax.gridlines()
    ax.add_feature(cfeature.LAND,facecolor=(0,1,0))
    ax.scatter(elon, elat, c='b',marker='*', s=275, transform=ccrs.PlateCarree())
    if len(w)>0:
        for i in range (len(w)):
            slon=w[i].stats.sac.stlo
            slat=w[i].stats.sac.stla
            ax.scatter(slon, slat, marker='o',c='r', s=75, transform=ccrs.PlateCarree())
            plt.text(slon - 2, slat - 2, w[i].stats.station,
                 horizontalalignment='right',
                 transform=ccrs.PlateCarree())
    else:
        for i in range(len(slat)):
            sln=slon[i]
            slt=slat[i]
            sta=stas[i]
            ax.scatter(sln, slt, marker='o',c='r', s=75, transform=ccrs.PlateCarree())
            plt.text(sln - 2, slt - 2, sta,
                 horizontalalignment='right',
                 transform=ccrs.PlateCarree())
    plt.show()