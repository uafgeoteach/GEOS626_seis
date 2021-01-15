import numpy as np
import spshell_love_config

def earthfun_love(r):
    """ Python adaptation of earthfun_love.m by Carl Tape
        Amanda McPherson, Jan 2021
        
        EARTHFUN_LOVE return a rho and mu value for a specified radius r

        called by stress_disp_love.py"""
    
    rspan = spshell_love_config.rspan
    crho = spshell_love_config.crho
    cmu = spshell_love_config.cmu
    mrho = spshell_love_config.mrho
    mmu = spshell_love_config.mmu
    cthick = spshell_love_config.cthick
    
    #ENTER CODE HERE
    
    return rho, mu