#!/usr/bin/env python

import sys,string,os
sys.path.insert(1, os.path.join(sys.path[0], '..')) ### import parent dir dependencies

import numpy as np
import numpy.polynomial.polynomial as poly
import matplotlib
import matplotlib.pyplot as plt

def disp(hgvfile):
    lab={}
    me=[]
    std=[]
    for line in open(hgvfile,'rU').xreadlines():
        data = line.rstrip()
        t = string.split(data,'\t')
        me.append(float(t[2]))
        std.append(float(t[1]))
        
    me=np.asarray(me)
    std=np.asarray(std)
    coefs = poly.polyfit(me, std, 2)
    ffit = poly.Polynomial(coefs)
    print len(ffit(me))
    np.savetxt("fitted.txt",ffit(me),delimiter="\t")
    fig1 = plt.figure()                                                                                           
    ax1 = fig1.add_subplot(111)                                                                                   
    ax1.scatter(me, std, facecolors='None')                                                                     
    ax1.plot(me, ffit(me))                                                                     
    plt.show()
    
disp("/Users/meenakshi/Documents/altanalyze-master/hgv_0.1.txt")
