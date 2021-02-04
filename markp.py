#!/usr/bin/env python
# coding: utf-8

# In[ ]:

# To use function: cid=fig.canvas.mpl_connect('button_press_event', markp)
# 
# This marks the 1/x value on a plot where you click the mouse.
# python version based on the markt.m writtin by Carl Tape
#
# INPUT ARGUMENTS:
# 
# left mouse button click: plot point
#
# EXAMPLE:
#
# fig = plt.figure()
# ax = fig.add_subplot(111)
# ax.set_xlim([0, 10])
# ax.set_ylim([0, 10]) 
# markt(fig)
#
# pt coordinates will print to screen

import matplotlib.pyplot as plt
import numpy as np
  
def markp(event):
    print('x=%f, y=%f' % (event.xdata, event.ydata))
    prd=round(1/event.xdata, 3)
    axe=event.inaxes
    axe.text(event.xdata, event.ydata, s = str(prd))
    #plt.plot(event.xdata, event.ydata, 'ro')
    plt.draw()