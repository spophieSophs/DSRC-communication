"""This script is the collision probabilistic model for DSRC technology."""
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fsolve
import math
import pickle
import os

def coll_prob(N1, N2, T):
    """This function calcualtes the direct collision probability by using analytical model"""
    """"
    Input:
        N1:     the number of neighbouring nodes around transmitting node, should be larger than 2
        N2:     the number of hidden terminals around receiving node, should be larger than 2
        T:      complete transmission time, in slots(13 microsec.)
    """
    # CW: CWmin contention window size, with user priority of 7, (aCWmin+1)/4-1 = 3, aCWmin is 15.
    CW = 3
    # Wm: average number of backoff slots preceding a transmission
    Wm = (CW+0)/2
    # tau: probability that a vehicle attempts to transmit a packet in an arbitary slot given that it has a packet in the queue
    tau = 1/(Wm +1)
    # sigma: duration of a backoff slot
    sigma = 1
    # E_T: average packet transmission time, assuming all packets have same length
    E_T = T
    # lamda: in packets per slot
    lamda = 13*(10**(-5))
    if N1 > 3:
        # this function solves the direct collision probability 
        def f(z):
            r_ = z[0]
            Pb = z[1]
            Pdc = z[2]
            K = z[3]
            E_S = z[4]
            f = np.zeros(5)
            # r_: model each station as an M/G/1/infinite buffer queue, r_ is the queue utilization i.e. r_ = lamda * E_S
            f[0] = lamda * E_S - r_
            # Pb: prob. that the channel is sensed busy when a new packet arrives, consider number of packets in one collision
            f[1] = (N1-1)*lamda*T*(1-Pdc/K*(K-1)) - Pb 
            # Pdc: prob. of direct collision i.e. Pdc = (1-(1-r_)*(1-Pb))*(1-(1-r_*tau)**(N-1))
            f[2] = (1-(1-r_)*(1-Pb))*(1-(1-r_*tau)**(N1-1)) - Pdc 
            # no: number of packets in one collisions
            f[3] = 2.0 + (1 - math.exp(-lamda*(N1-2)*T)) *lamda*(N1-2)*T- K
            # E_Y: average interruption period per slot, for simplicity, assume that every backoff slot can be interrupted at most once
            E_Y = (1- (1- r_ * tau)**(N1-1)) * T
            # E_U: average contention window size
            E_U = Wm 
            # E_B: average total backoff duration
            E_B = (sigma + E_Y) * E_U
            # E_Tres: average residual lifetime of an ongoing transmission Tres
            E_Tres = 1/2*(lamda*T-1+math.exp(-lamda*T))/((1-1/2*math.exp(-lamda*T))*lamda)
            # E_A: average access delay
            E_A = (1- r_)* Pb * (E_B + E_Tres) + r_ * E_B
            # E_S: average service time, the combination of channel access delay and transmission time
            f[4] = E_A + E_T - E_S            
            return f
        z = fsolve(f, [0.01, 0.4, 0.2, 2, 50])
        # DC_prob: direct collision probability
        DC_prob = z[2]
    else:
        DC_prob = 1
    if N2 > 3:
        # this function solves the direct collision probability 
        def f(z):
            r_ = z[0]
            Pb = z[1]
            Pdc = z[2]
            K = z[3]
            E_S = z[4]
            f = np.zeros(5)
            # r_: model each station as an M/G/1/infinite buffer queue, r_ is the queue utilization i.e. r_ = lamda * E_S
            f[0] = lamda * E_S - r_
            # Pb: prob. that the channel is sensed busy when a new packet arrives, consider number of packets in one collision
            f[1] = (N2-1)*lamda*T*(1-Pdc/K*(K-1)) - Pb 
            # Pdc: prob. of direct collision i.e. Pdc = (1-(1-r_)*(1-Pb))*(1-(1-r_*tau)**(N-1))
            f[2] = (1-(1-r_)*(1-Pb))*(1-(1-r_*tau)**(N2-1)) - Pdc 
            # no: number of packets in one collisions
            f[3] = 2.0 + (1 - math.exp(-lamda*(N2-2)*T)) *lamda*(N2-2)*T- K
            # E_Y: average interruption period per slot, for simplicity, assume that every backoff slot can be interrupted at most once
            E_Y = (1- (1- r_ * tau)**(N2-1)) * T
            # E_U: average contention window size
            E_U = Wm 
            # E_B: average total backoff duration
            E_B = (sigma + E_Y) * E_U
            # E_Tres: average residual lifetime of an ongoing transmission Tres
            E_Tres = 1/2*(lamda*T-1+math.exp(-lamda*T))/((1-1/2*math.exp(-lamda*T))*lamda)
            # E_A: average access delay
            E_A = (1- r_)* Pb * (E_B + E_Tres) + r_ * E_B
            # E_S: average service time, the combination of channel access delay and transmission time
            f[4] = E_A + E_T - E_S            
            return f
        z = fsolve(f, [0.01, 0.4, 0.2, 2, 50])
        # DC_prob: direct collision probability
        HDC_prob = z[2]
    else: 
        HDC_prob = 1
    if N2 > 1: 
        # calculate the packets number in collision of hidden terminal nodes
        Kcol = 2 + (lamda*(N2-1)*T) *(1-math.exp(-lamda*(N2-2)*T))
        # the first event that needs to be fulfilled of a successful non-hidden terminal collision
        ht1 = 1 - (N2-1)* lamda * T *(1-HDC_prob/Kcol * (Kcol-1))
        # the second event that needs to be fulfilled of a successful non-hidden terminal collision
        ht2 = math.exp(-lamda*N2*(T+2))
        # the packet reception probability calculation
        PR_prob =(1- DC_prob)*ht1*ht2 
        # if the hidden terminal number is less than 1, no need for hidden terminal check.
    else:
        PR_prob = 1 - DC_prob
    return PR_prob
