# -*- coding: utf-8 -*-
"""
Created on Tue Apr 12 00:18:47 2016

@author: Sergey
"""

from astropy.io import fits
import numpy as NP
from scipy import ndimage
from astropy import coordinates
from astropy import constants
from BadaryRAO import BadaryRAO

class SrhFitsFile():
    def __init__(self, name, sizeOfUv):
        self.omegaEarth = coordinates.earth.OMEGA_EARTH.to_value()
        self.isOpen = False;
        self.calibIndex = 10;
        self.frequencyChannel = 8;
        self.centerP = 0.;
        self.deltaP = 0.005;
        self.centerQ = 0.;
        self.deltaQ = 0.005;
        self.centerH = 0.;
        self.deltaH = 4.9/3600.*NP.pi;
        self.centerD = 0.;
        self.deltaD = 4.9/3600.*NP.pi;
        self.ewLcpPhaseCorrection = NP.zeros(32)
        self.ewRcpPhaseCorrection = NP.zeros(32)
        self.sLcpPhaseCorrection = NP.zeros(16)
        self.sRcpPhaseCorrection = NP.zeros(16)
        self.ewPhaMatrix = [ \
            [1, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1]];
        self.ewAmpMatrix = [ \
            [1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1]];
        self.sPhaMatrix = [ \
            [1, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1]];
        self.sAmpMatrix = [ \
            [1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1]];
        self.ewPhaseClosureMatrix = [ \
            [ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [ 0,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [ 0, 0,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [ 0, 0, 0,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [ 0, 0, 0, 0,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [ 0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [ 0, 0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [ 0, 0, 0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [ 0, 0, 0, 0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [ 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0], \
            [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0, 0, 0, 0], \
            [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0, 0, 0], \
            [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0, 0], \
            [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0], \
            [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 1, 0, 0, 0], \
            [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 1, 0, 0], \
            [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 1, 0], \
            [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 1]]

        self.sPhaseClosureMatrix = [ \
            [ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [ 0,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [ 0, 0,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [ 0, 0, 0,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [ 0, 0, 0, 0,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [ 0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [ 0, 0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0], \
            [ 0, 0, 0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0, 0, 0, 0], \
            [ 0, 0, 0, 0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0, 0, 0], \
            [ 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0, 0], \
            [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0], \
            [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 1, 0, 0, 0], \
            [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 1, 0, 0], \
            [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 1, 0], \
            [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 1]];
        self.sizeOfUv = sizeOfUv
        self.open(name)
                                    
    def open(self,name):
        try:
            self.hduList = fits.open(name);
            self.isOpen = True;
            fitsDate = self.hduList[0].header['DATE-OBS'].split('/')
            fitsTime = self.hduList[0].header['TIME-OBS'].split('.')
            fitsTime = fitsTime[0].split(':')
            self.dateObs = fitsDate[0] + '-' + fitsDate[1] + '-' + fitsDate[2] + 'T' \
                            + fitsTime[0] + ':' + fitsTime[1] + ':' + fitsTime[2]
            self.antennaNumbers = self.hduList[2].data['ANTENNA'];
            self.antennaNumbers = NP.reshape(self.antennaNumbers,self.antennaNumbers.size);
            self.antennaA = self.hduList[2].data['ANTA'];
            self.antennaA = NP.reshape(self.antennaA,self.antennaA.size);
            self.antennaB = self.hduList[2].data['ANTB'];
            self.antennaB = NP.reshape(self.antennaB,self.antennaB.size);
            self.uvLcp = NP.zeros((self.sizeOfUv,self.sizeOfUv),dtype=complex);
            self.uvRcp = NP.zeros((self.sizeOfUv,self.sizeOfUv),dtype=complex);
            self.indexA = (self.antennaA - 1)*2;
            self.indexB = (self.antennaB - 1 - 128)*2;
            self.freqList = self.hduList[1].data['FREQUENCY'];
            self.freqListLength = self.freqList.size;
            self.dataLength = self.hduList[1].data['TIME'].size // self.freqListLength;
#            self.freqTime = NP.reshape(self.hduList[1].data['TIME'],(self.freqListLength,self.dataLength));
            self.freqTime = self.hduList[1].data['TIME']
            self.visListLength = self.hduList[1].data['VIS_LCP'].size // self.freqListLength // self.dataLength;
            self.visLcp = NP.reshape(self.hduList[1].data['VIS_LCP'],(self.freqListLength,self.dataLength,self.visListLength));
            self.visRcp = NP.reshape(self.hduList[1].data['VIS_RCP'],(self.freqListLength,self.dataLength,self.visListLength));
            self.ampLcp = NP.reshape(self.hduList[1].data['AMP_LCP'],(self.freqListLength,self.dataLength,self.antennaNumbers.size));
            self.ampRcp = NP.reshape(self.hduList[1].data['AMP_RCP'],(self.freqListLength,self.dataLength,self.antennaNumbers.size));
            self.visLcp.real = NP.sin(NP.pi/2*self.visLcp.real);
            self.visRcp.real = NP.sin(NP.pi/2*self.visRcp.real);
            self.visLcp.imag = NP.sin(NP.pi/2*self.visLcp.imag);
            self.visRcp.imag = NP.sin(NP.pi/2*self.visRcp.imag);

            self.RAO = BadaryRAO(self.dateObs.split('T')[0])

        except FileNotFoundError:
            print('File %s  not found'%name);
    
    def append(self,name):
        try:
            hduList = fits.open(name);
#            freqTime = NP.reshape(self.hduList[1].data['TIME'],(self.freqListLength,self.dataLength));
            freqTime = hduList[1].data['TIME']
            dataLength = hduList[1].data['TIME'].size // self.freqListLength;
            visLcp = NP.reshape(hduList[1].data['VIS_LCP'],(self.freqListLength,dataLength,self.visListLength));
            visRcp = NP.reshape(hduList[1].data['VIS_RCP'],(self.freqListLength,dataLength,self.visListLength));
            ampLcp = NP.reshape(hduList[1].data['AMP_LCP'],(self.freqListLength,dataLength,self.antennaNumbers.size));
            ampRcp = NP.reshape(hduList[1].data['AMP_RCP'],(self.freqListLength,dataLength,self.antennaNumbers.size));
            visLcp.real = NP.sin(NP.pi/2*visLcp.real);
            visRcp.real = NP.sin(NP.pi/2*visRcp.real);
            visLcp.imag = NP.sin(NP.pi/2*visLcp.imag);
            visRcp.imag = NP.sin(NP.pi/2*visRcp.imag);

            self.freqTime = NP.concatenate((self.freqTime, freqTime), axis = 1)
            self.visLcp = NP.concatenate((self.visLcp, visLcp), axis = 1)
            self.visRcp = NP.concatenate((self.visRcp, visRcp), axis = 1)
            self.ampLcp = NP.concatenate((self.ampLcp, ampLcp), axis = 1)
            self.ampRcp = NP.concatenate((self.ampRcp, ampRcp), axis = 1)
            self.dataLength += dataLength
            hduList.close()

        except FileNotFoundError:
            print('File %s  not found'%name);

    def getHourAngle(self, scan):
        self.hAngle = self.omegaEarth * (self.freqTime[self.frequencyChannel, scan] - self.RAO.culmination)
        return self.hAngle
        
    def setSizeOfUv(self, sizeOfUv):
        self.sizeOfUv = sizeOfUv
        self.uvLcp = NP.zeros((self.sizeOfUv,self.sizeOfUv),dtype=complex);
        self.uvRcp = NP.zeros((self.sizeOfUv,self.sizeOfUv),dtype=complex);
        
    def getDeclination(self):
        return self.RAO.declination

    def getPQScale(self, size, FOV):
        self.cosP = NP.sin(self.hAngle) * NP.cos(self.RAO.declination)
        self.cosQ = NP.cos(self.hAngle) * NP.cos(self.RAO.declination) * NP.sin(self.RAO.observatory.lat) - NP.sin(self.RAO.declination) * NP.cos(self.RAO.observatory.lat)
        FOV_p = 2.*(constants.c / (self.freqList[self.frequencyChannel]*1e6)) / (self.RAO.base*NP.sqrt(1. - self.cosP**2.));
        FOV_q = 2.*(constants.c / (self.freqList[self.frequencyChannel]*1e6)) / (self.RAO.base*NP.sqrt(1. - self.cosQ**2.));
        
        return [int(size*FOV/FOV_p.to_value()), int(size*FOV/FOV_q.to_value())]
        
    def getPQ2HDMatrix(self):
        gP =  NP.arctan(NP.tan(self.hAngle)*NP.sin(self.RAO.declination));
        gQ =  NP.arctan(-(NP.sin(self.RAO.declination) / NP.tan(self.hAngle) + NP.cos(self.RAO.declination) / (NP.sin(self.hAngle)*NP.tan(self.RAO.observatory.lat))));
        
        if self.hAngle > 0:
            gQ = NP.pi + gQ;
        g = gP - gQ;
          
        pqMatrix = NP.zeros((3,3))
        pqMatrix[0, 0] =  NP.cos(gP) - NP.cos(g)*NP.cos(gQ)
        pqMatrix[0, 1] = -NP.cos(g)*NP.cos(gP) + NP.cos(gQ)
        pqMatrix[1, 0] =  NP.sin(gP) - NP.cos(g)*NP.sin(gQ)
        pqMatrix[1, 1] = -NP.cos(g)*NP.sin(gP) + NP.sin(gQ)
        pqMatrix /= NP.sin(g)**2.
        pqMatrix[2, 2] = 1.
        return pqMatrix
        
    def close(self):
        self.hduList.close();
    
    def updateAntennaAmplitude(self):
        self.ewAmpLcp = NP.abs(self.visLcp[self.frequencyChannel,self.calibIndex,512+15:512+15+31]);
        self.ewAmpRcp = NP.abs(self.visRcp[self.frequencyChannel,self.calibIndex,512+15:512+15+31]);
        self.sAmpLcp  = NP.abs(self.visLcp[self.frequencyChannel,self.calibIndex,512:512+15]);
        self.sAmpRcp  = NP.abs(self.visRcp[self.frequencyChannel,self.calibIndex,512:512+15]);

#        self.ewAntAmpLcp, c, d, e = NP.linalg.lstsq(self.ewAmpMatrix,NP.log(self.ewAmpLcp))
#        self.ewAntAmpRcp, c, d, e = NP.linalg.lstsq(self.ewAmpMatrix,NP.log(self.ewAmpRcp))
        self.ewAntAmpLcp0, c, d, e = NP.linalg.lstsq(self.sAmpMatrix,NP.log(self.ewAmpLcp[0:15]))
        self.ewAntAmpRcp0, c, d, e = NP.linalg.lstsq(self.sAmpMatrix,NP.log(self.ewAmpRcp[0:15]))
        self.ewAntAmpLcp1, c, d, e = NP.linalg.lstsq(self.sAmpMatrix,NP.log(self.ewAmpLcp[16:31]))
        self.ewAntAmpRcp1, c, d, e = NP.linalg.lstsq(self.sAmpMatrix,NP.log(self.ewAmpRcp[16:31]))
        
        self.ewAntAmpLcp = NP.zeros(33)
        self.ewAntAmpRcp = NP.zeros(33)
        self.ewAntAmpLcp[1:16] = self.ewAntAmpLcp0[1:16]
        self.ewAntAmpLcp[17:32] = self.ewAntAmpLcp1[1:16]
        self.ewAntAmpRcp[1:16] = self.ewAntAmpRcp0[1:16]
        self.ewAntAmpRcp[17:32] = self.ewAntAmpRcp1[1:16]

        self.sAntAmpLcp, c, d, e = NP.linalg.lstsq(self.sAmpMatrix,NP.log(self.sAmpLcp))
        self.sAntAmpRcp, c, d, e = NP.linalg.lstsq(self.sAmpMatrix,NP.log(self.sAmpRcp))

        self.ewAntAmpLcp = NP.exp(self.ewAntAmpLcp)
        self.ewAntAmpRcp = NP.exp(self.ewAntAmpRcp)
        self.sAntAmpLcp  = NP.exp(self.sAntAmpLcp)
        self.sAntAmpRcp  = NP.exp(self.sAntAmpRcp)

        self.ewAntAmpLcp = self.ewAntAmpLcp / NP.max(self.ewAntAmpLcp);
        self.ewAntAmpRcp = self.ewAntAmpRcp / NP.max(self.ewAntAmpRcp);
        self.sAntAmpLcp  = self.sAntAmpLcp / NP.max(self.sAntAmpLcp);
        self.sAntAmpRcp  = self.sAntAmpRcp / NP.max(self.sAntAmpRcp);

    def updateAntennaPhase(self):
        self.ewPhaLcp = NP.angle(self.visLcp[self.frequencyChannel,self.calibIndex,512+15:512+15+31]);
        self.ewPhaRcp = NP.angle(self.visRcp[self.frequencyChannel,self.calibIndex,512+15:512+15+31]);
        self.sPhaLcp  = NP.angle(self.visLcp[self.frequencyChannel,self.calibIndex,512:512+15]);
        self.sPhaRcp  = NP.angle(self.visRcp[self.frequencyChannel,self.calibIndex,512:512+15]);

        self.ewAntPhaLcp, c, d, e = NP.linalg.lstsq(self.ewPhaMatrix,self.ewPhaLcp);
        self.sAntPhaLcp, c, d, e = NP.linalg.lstsq(self.sPhaMatrix,self.sPhaLcp);
        self.ewAntPhaRcp, c, d, e = NP.linalg.lstsq(self.ewPhaMatrix,self.ewPhaRcp);
        self.sAntPhaRcp, c, d, e = NP.linalg.lstsq(self.sPhaMatrix,self.sPhaRcp);
        
    def setCalibIndex(self, calibIndex):
        self.calibIndex = calibIndex;
        self.updateAntennaPhase()
        self.updateAntennaAmplitude()

    def setFrequencyChannel(self, channel):
        self.frequencyChannel = channel
        self.updateAntennaPhase()
        self.updateAntennaAmplitude()
        
    def changeEastWestPhase(self, newLcpPhaseCorrection, newRcpPhaseCorrection):
        self.ewLcpPhaseCorrection[:] = newLcpPhaseCorrection[:]
        self.ewRcpPhaseCorrection[:] = newRcpPhaseCorrection[:]
        
    def changeSouthPhase(self, newLcpPhaseCorrection, newRcpPhaseCorrection):
        self.sLcpPhaseCorrection[:] = newLcpPhaseCorrection[:]
        self.sRcpPhaseCorrection[:] = newRcpPhaseCorrection[:]

    def vis2uv(self,  scan, phaseCorrect = False, amplitudeCorrect = False, PSF = False):
        self.uvLcp[:,:] = complex(0,0)
        self.uvRcp[:,:] = complex(0,0)
        O = self.sizeOfUv//2 + 1
        if PSF == False:
            for jj in range(16):
                for ii in range(32):
                  self.uvLcp[O + jj*2, O + (ii - 16)*2] = self.visLcp[self.frequencyChannel, scan, jj*32 + ii]
                  self.uvRcp[O + jj*2, O + (ii - 16)*2] = self.visRcp[self.frequencyChannel, scan, jj*32 + ii]
                  self.uvLcp[O - 2 - jj*2, O + (15 - ii)*2] = NP.conj(self.uvLcp[O + jj*2, O + (ii - 16)*2])
                  self.uvRcp[O - 2 - jj*2, O + (15 - ii)*2] = NP.conj(self.uvRcp[O + jj*2, O + (ii - 16)*2])
            if (phaseCorrect):
                for jj in range(16):
                    for ii in range(32):
                        self.uvLcp[O + jj*2, O + (ii - 16)*2] *= NP.exp(1j*((self.ewAntPhaLcp[ii + 1] + self.ewLcpPhaseCorrection[ii]) - (self.sAntPhaLcp[jj + 1]) + self.sLcpPhaseCorrection[jj]))
                        self.uvRcp[O + jj*2, O + (ii - 16)*2] *= NP.exp(1j*((self.ewAntPhaRcp[ii + 1] + self.ewRcpPhaseCorrection[ii]) - (self.sAntPhaRcp[jj + 1]) + self.sRcpPhaseCorrection[jj]))
                        self.uvLcp[O - 2 - jj*2, O + (15 - ii)*2] = NP.conj(self.uvLcp[O + jj*2, O + (ii - 16)*2])
                        self.uvRcp[O - 2 - jj*2, O + (15 - ii)*2] = NP.conj(self.uvRcp[O + jj*2, O + (ii - 16)*2])
            if (amplitudeCorrect):
                for jj in range(16):
                    for ii in range(32):
                        self.uvLcp[O + jj*2, O + (ii - 16)*2] /= (self.ewAntAmpLcp[ii + 1] * self.sAntAmpLcp[jj + 1])
                        self.uvRcp[O + jj*2, O + (ii - 16)*2] /= (self.ewAntAmpRcp[ii + 1] * self.sAntAmpRcp[jj + 1])
                        self.uvLcp[O - 2 - jj*2, O + (15 - ii)*2] = NP.conj(self.uvLcp[O + jj*2, O + (ii - 16)*2])
                        self.uvRcp[O - 2 - jj*2, O + (15 - ii)*2] = NP.conj(self.uvRcp[O + jj*2, O + (ii - 16)*2])
        else:
            for jj in range(16):
                for ii in range(32):
                  self.uvLcp[O + jj*2, O + (ii - 16)*2] = .01 + 0j
                  self.uvRcp[O + jj*2, O + (ii - 16)*2] = .01 + 0j
                  self.uvLcp[O - 2 - jj*2, O + (15 - ii)*2] = NP.conj(self.uvLcp[O + jj*2, O + (ii - 16)*2])
                  self.uvRcp[O - 2 - jj*2, O + (15 - ii)*2] = NP.conj(self.uvRcp[O + jj*2, O + (ii - 16)*2])
            if (phaseCorrect):
                for jj in range(16):
                    for ii in range(32):
                        self.uvRcp[O + jj*2, O + (ii - 16)*2] *= NP.exp(1j*((-self.ewAntPhaRcp[ii + 1] - self.ewRcpPhaseCorrection[ii]) + (self.sAntPhaRcp[jj + 1]) + self.sRcpPhaseCorrection[jj]))
                        self.uvRcp[O - 2 - jj*2, O + (15 - ii)*2] = NP.conj(self.uvRcp[O + jj*2, O + (ii - 16)*2])
            if (amplitudeCorrect):
                for jj in range(16):
                    for ii in range(32):
                        self.uvRcp[O + jj*2, O + (ii - 16)*2] *= (self.ewAntAmpRcp[ii + 1] * self.sAntAmpRcp[jj + 1])
                        self.uvRcp[O - 2 - jj*2, O + (15 - ii)*2] = NP.conj(self.uvRcp[O + jj*2, O + (ii - 16)*2])

    def uv2lmImage(self):
        self.lcp = NP.fft.fft2(NP.roll(NP.roll(self.uvLcp,self.sizeOfUv//2,0),self.sizeOfUv//2,1));
        self.lcp = NP.roll(NP.roll(self.lcp,self.sizeOfUv//2,0),self.sizeOfUv//2,1);
        self.rcp = NP.fft.fft2(NP.roll(NP.roll(self.uvRcp,self.sizeOfUv//2,0),self.sizeOfUv//2,1));
        self.rcp = NP.roll(NP.roll(self.rcp,self.sizeOfUv//2,0),self.sizeOfUv//2,1);
        
    def ij2hd(self, ij):
        return(((ij[0] - self.centerH)*self.deltaH,(ij[1] - self.centerD)*self.deltaD));

    def hd2pq(self, hd):
        return((hd[0],hd[1]));

    def pq2kl(self, pq):
        return((pq[0]/self.deltaP + self.centerP,pq[1]/self.deltaQ + self.centerQ));

    def ij2kl(self, ij):
        return(self.pq2kl(self.hd2pq(self.ij2hd(ij))));
        
    def lmImage2hdImage(self):
        return(ndimage.geometric_transform(self.lcp.real,self.ij2kl));
    
    def getDateObs(self):
        return self.dateObs if self.isOpen else ''

    def getDataLength(self):
        return self.dataLength if self.isOpen else 0

    def getTimesObs(self):
        return self.hduList[1].data['TIME'] if self.isOpen else 0

    def getAntennaA(self):
        return self.antennaA if self.isOpen else 0
        
    def getAntennaB(self):
        return self.antennaB if self.isOpen else 0

    def getVisLcp(self):
        return self.visLcp if self.isOpen else 0
        
    def getVisRcp(self):
        return self.visRcp if self.isOpen else 0
    
    def getFrequencyList(self):
        return self.freqList
    
    def phaseClosure(self, antA, antB, antC):
        pC = NP.zeros(self.dataLength)
        for i in range(self.dataLength):
            pC[i] = NP.angle(self.visLcp[5,i,antA]) + NP.angle(self.visLcp[5,i,antB]) - NP.angle(self.visLcp[5,i,antC])
        return pC

    def phaseAnt(self, antA):
        pC = NP.zeros(self.dataLength)
        for i in range(self.dataLength):
            pC[i] = NP.angle(self.visLcp[5,i,antA]) 
        return pC

    def magnitAnt(self, antA):
        pC = NP.zeros(self.dataLength)
        for i in range(self.dataLength):
            pC[i] = NP.abs(self.visLcp[5,i,antA]) 
        return pC

    def phaseClosureSouthVector(self, visEw, scan, freq):
        pC = NP.zeros(16)
        for i in range(15):
            pC[i + 1] = NP.angle(self.visLcp[freq,scan,512 + 15 + visEw]) + NP.angle(self.visLcp[freq,scan,32*i + visEw + 1]) - NP.angle(self.visLcp[freq,scan,32*i + visEw])
        return pC

    def phaseClosureEastWestVector(self, visS, scan, freq):
        pC = NP.zeros(32)
        pC[0] = 1.
        for i in range(31):
            pC[i + 1] = - NP.angle(self.visLcp[freq,scan,16*visS + i]) + NP.angle(self.visLcp[freq,scan,512 + visS]) - NP.angle(self.visLcp[freq,scan,16*visS + i + 1]) 
        return pC
    
    def solveSouthPhaseClosure(self, pcVector):
        phases, c, d, e = NP.linalg.lstsq(self.sPhaseClosureMatrix,pcVector);
        return phases

    def solveEastWestPhaseClosure(self, pcVector):
        phases, c, d, e = NP.linalg.lstsq(self.ewPhaseClosureMatrix,pcVector);
        return phases
    