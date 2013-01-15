#################################################################
# @Program: buildTFBSDB.py                                      #
# @Version: 1                                                   #
# @Author: Christopher L Plaisier, PhD                          #
# @Sponsored by:                                                #
# Nitin Baliga, ISB                                             #
# Institute for Systems Biology                                 #
# 401 Terry Ave North                                           #
# Seattle, Washington  98109-5234                               #
# (216) 732-2139                                                #
# @Also Sponsored by:                                           #
# Luxembourg Systems Biology Grant                              #
# American Cancer Society Postdoctoral Fellowship               #
#                                                               #
# If this program is used in your analysis please mention who   #
# built it. Thanks. :-)                                         #
#                                                               #
# Copyrighted by Chris Plaisier  1/15/2013                      #
#################################################################

###############
### IMPORTS ###
###############
from pssm import pssm
import cPickle, gzip, os, sys, re, os, math, shutil
from copy import deepcopy
from subprocess import *
from random import sample
from multiprocessing import Pool, cpu_count, Manager

###
# 
###



prev = ['chr0',0,0]
chrs = {}
inFile = open('combined.fps','r')
while 1:
    inLine = inFile.readline()
    
inFile.close()
