
import numpy as np
import matplotlib.pyplot as plt
from sph2cart import sph2cart
from itertools import product, combinations

def globefun3(R,lat,lon,bool_point,lc, fig,ax):
    # Python adaptation of globefun3.m
    # Python coding done by Nealey Sims
    #
    # GLOBEFUN3 plots a single point on the sphere with a reference globe
    #
    # INPUT
    #   R           radius of sphere
    #   lat         latitude (deg)  (theta = 90 - lat)
    #   lon         longitude (deg) (lon = phi)
    #   bool_point  plot the specified point
    #   lc          line style/color
    #
    #
    
    deg = 180/np.pi
    lw = 1.0
    th = (90-lat)/deg
    ph = lon/deg
    
    # (1) xy-plane circle (equator)
    # (2) yz-plane
    # (3) xz-plane (prime meridian)
    gpts = 100
    tt = np.linspace(0, 2*np.pi, gpts)
    
    
    ax.set_box_aspect((4,4,3.4))  ## Equal axis
    ax.plot(R*np.cos(tt), R*np.sin(tt), np.zeros(gpts), lc, linewidth=lw)
    ax.plot(np.zeros(gpts), R*np.cos(tt), R*np.sin(tt), lc, linewidth=lw)
    ax.plot(R*np.cos(tt), np.zeros(gpts), R*np.sin(tt), lc, linewidth=lw)
    
    fac = 1.15; lw2 = 0.5
    RR = np.array([-R,R])
    ax.plot((fac*RR),[0, 0],[0, 0],'k', linewidth=lw2)
    ax.plot([0, 0],(fac*RR),[0, 0],'k', linewidth=lw2)
    ax.plot([0, 0],[0, 0],(fac*RR),'k', linewidth=lw2)
    
    # draw cube
    r = (fac*RR)
    for s, e in combinations(np.array(list(product(r, r, r))), 2):
        if np.sum(np.abs(s-e)) == r[1]-r[0]:
            ax.plot3D(*zip(s, e), color="k", linewidth=lw2)
    
    
    if bool_point == 1:
        # plot the specified point P(r,th,phi)
        msize = 20; lw3 = lw+0.5
        # [X,Y,Z] = sph2cart(az,el,r)
        [X,Y,Z] = sph2cart(ph,np.pi/2-th,R)  #Note: need sph2cart function in folder
        
        [X2,Y2,Z2] = sph2cart(ph,0,R)
  
        ax.plot([X],[Y],[Z], color=lc, marker='o')
        ax.plot([0, X],[0, Y],[0, Z],lc, linewidth= lw3)
        
        leps = 0.1;
        if ( lat >= -90+leps) and (lat <= 90-leps):
            xzplane = np.array([R*np.cos(tt) , np.zeros(gpts) , R*np.sin(tt)])
            rotmat = np.array([[np.cos(ph), -np.sin(ph), 0] , [np.sin(ph), np.cos(ph), 0], [0, 0, 1]])
            rotxz = np.dot(rotmat,xzplane)
            ax.plot(rotxz[0,:], rotxz[1,:], rotxz[2,:],lc,linewidth=lw);
            aa = lat/deg;
            ax.plot(R*np.cos(aa)*np.cos(tt), R*np.cos(aa)*np.sin(tt), R*np.sin(aa)*np.ones(gpts),lc,linewidth =lw);
            
            # plots some extra lines for reference
            ax.plot([0, X2],[0, Y2],[0, 0],lc, linewidth= lw2);
            ax.plot([X, X],[Y, Y],[0, Z],lc, linewidth= lw2);
            ax.plot([0, X],[0, Y],[Z, Z],lc, linewidth= lw2);
        #plt.title('\u03C6 = %.2f, \u03B8 = %.2f' % (ph*deg, th*deg))
        ax.view_init(10, 10)
    ax.grid(False)
    ax.xaxis.set_pane_color((0.0, .0, 0.0, 0.0))
    ax.yaxis.set_pane_color((0.0, .0, 0.0, 0.0))
    ax.zaxis.set_pane_color((0.0, .0, 0.0, 0.0))
    #plt.show()