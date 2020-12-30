import numpy as np
import spshell_config

def earthfun(r):
    """ Python adaptation of earthfun.m by Carl Tape
        Amanda McPherson, Dec 2020
        
        EARTHFUN return a rho and mu value for a specified radius r
        Note that rho and mu are NOT globally defined within this function;
        therefore what happens to rho and mu inside here will NOT affect the
        global values of rho and mu.

        called by stress_disp_tor.m"""
    
    #global rspan, imod
    imod = spshell_config.imod
    rspan = spshell_config.rspan
    
    if imod == 1:
        # linear model
        # enter your code here
        cmbr = rspan[0]    # b
        earthr = rspan[1]  # a
        
        raise Exception('earthfun.py imod=1 not yet implemented (comment this out when you have implemented it)')
        #spshell_config.rho = rho  #uncomment when you've implemented your code
        #spshell_config.mu = mu  #uncomment when you've implemented your code
        
    elif imod == 2:
        # cubic model
        # enter your code here
        
        raise Exception('earthfun.py imod=2 not yet implemented (comment this out when you have implemented it)')
        #spshell_config.rho = rho  #uncomment when you've implemented your code
        #spshell_config.mu = mu  #uncomment when you've implemented your code
        
    elif imod != 1 or imod != 2:
        raise Exception('invalid imod (=1,2)')
        
    return rho, mu