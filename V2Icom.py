import os, sys, time, subprocess, pickle  
sys.path.append(os.path.join('c://','Program Files','sumo-0.30.0','tools'))
import traci
import math
from DeLOSNLOS import * 
from collisionprobability import coll_prob
import time, random
from rtree import index
import pandas as pd
import itertools
"""Update version, includes rec_pow, positions, collision probability, and etc. Excluding the number of neighbours.
   NOTE-1: the position of vehicles has been changed based on veins rule-- add margin, and y flipflop.
   NOTE-2: the rec_pow function has been modified, thus comment out the antenna gain part.
   NOTE-3: the dist_thre has been modified, based on two ray interference model.
   NOTE-4: move the dist_thre lines to up, in order to speed up a bit.
   NOTE-5: to speed up the model, add neighbours to each vehicle before calculate the rec.pow.
   NOTE-6: the x-y coordinates should be the origin version, but change them to fit veins positions for output.
   NOTE-7: updated building geometry file, the previous version missed a lot of building geometries.
   NOTE-8: to do sensitivity analysis, changes the vehicle density. 5 different simulation time, each time with 10 seconds.
   NOTE-9: fix the niehgbours and hidden terminals part, add limitation of receiving power level.
   NOTE-10:for V2I communication. A permanent receiver: traffic light location."""


# input file: region of interest, use the xml2geojson file 
infile = r'path/to/building/geometry/geojsonfile'
# ROI: region of interest, x-y coordinates, left bottom and right up.
ROI = [(4894.00, 6732.00), (7395.00, 8175.00)]

# read building geometry file
with open(infile, 'r') as f:
    js = json.load(f)
# initiate the rtree index
idx = index.Index()
# insert index of building bounds
for pos, poly in enumerate(js['features']):
    idx.insert(pos, shape(poly['geometry']).bounds)

# antenna gain, in dB i
Gt = 1
Gr = 1
# antenna height of private vehicles, in m
txHeight = 1.895
rxHeight = 1.895
# antenna height of bus, in m
busHeight = 3.15
# antenna height of intersection, in m
intHeight = 5
# transmission power, in mW
Pt = 20
# receiving power threshold, in dBm
Rec_th = -89
# transmission distance threshold, in m. This should depend on if the power exceeds the threshold, -89dBm, used two ray interference path loss model
dist_th = 977
# packet transmission time, in slot
T_pk = 40

# port of traci connection to python
PORT = 8813
# sumoBinary = "path/to/sumo-gui.exe"
sumoBinary = "path/to/sumo.exe"
sumoConfig = "path/to/V2Isumocfg" # change the packet sending frequency at V2I.sumocfg and the step range in this code
sumoProcess = subprocess.Popen([sumoBinary, '-c', sumoConfig, "--remote-port", str(PORT)], stdout= sys.stdout, stderr=sys.stderr)
# setup traci
traci.init(PORT)


# initialize an empty list
busls = []
# initialize empty lists to contain number of vehicles in each step, number of buses in each step
Num_veh = []
Num_bus = []
# intersection location
intloc = [(6940, 6902), (7074, 7257), (6971, 7602)]

for step in range(5000):
    # warm up sumo simulation
    traci.simulationStep()
    if step < 4900:
        continue
    if step == 4901:
    # record the time of simulation 
        start_time = time.time()
    # initialize an empty dictionary containing the vehicle ids and positions
    default_xy = {}

    # get all vehicle id list
    vehIDlist = traci.vehicle.getIDList()
    # sumo simulation time
    curr_time = traci.simulation.getCurrentTime()
    for vehid in vehIDlist:
        # get vehicle x-y position, and write in a dictionary, add Veins' werid margin value
        x, y = traci.vehicle.getPosition(vehid)
        x = round(x, 2)
        y = round(y, 2)
        # check if the vehicle's position is in the region of interest area (use true values without margin)
        if meetsROIareaxml((x, y), ROI):
            default_xy.update({vehid:(x, y)})
    xypositions = list(default_xy.values()) 

    # number of vehicles in the ROI of each step
    Num_veh.append(len(xypositions))

    # check communication between buses and intersections
    for intersec in intloc:
        intNeigh = []
        # initialize an empty dictionary containing Rec.power between buses and infrastructure.
        busdef = {}
        recpos = intersec
        # get the neighbours of intersection, including buses
        for veh, vehloc in default_xy.items():
            vehint_dist = caldistance(vehloc, recpos)
            if vehint_dist > dist_th:
                continue
            walld_vehint, wallnum_vehint = LOSNLOSxml(js, idx, vehloc, recpos)
            if wallnum_vehint > 0:
                Rec_power_vehint = NonLineofSight(Pt, txHeight, intHeight, Gt, Gr, vehint_dist, walld_vehint, wallnum_vehint)
            else:
                Rec_power_vehint = LineofSight(Pt, txHeight, intHeight, Gt, Gr, vehint_dist) 
            if Rec_power_vehint > Rec_th:
                intNeigh.append(veh)
            if traci.vehicle.getVehicleClass(veh) == 'bus':
                busdef.update({veh: (vehloc, vehint_dist, Rec_power_vehint, walld_vehint, wallnum_vehint)})
        for veh, data in busdef.items():  
            busNeigh = []
            # Get the total length of intersection and the number of crossing walls of the building, x-y coordinates
            bus = veh
            busloc = data[0]
            TxRx_dist = data[1]
            # get the receiving power level between bus and infrasturcture
            Rec_power = data[2]
            walld = data[3]
            wallnum = data[4]
            if Rec_power > Rec_th:
                # vehtup includes the communicating pair that has a higher receiving power lvl than the threshold
                for veh, vehloc in default_xy.items():
                    if veh == bus:
                        continue
                    busprv_dist = caldistance(busloc, vehloc)
                    # ignore comm_pair with distance>977, in this way also speed up model and the results would be correct.
                    if busprv_dist > dist_th:
                        continue
                    walld_busprv, wallnum_busprv = LOSNLOSxml(js, idx, busloc, vehloc)
                    if wallnum_busprv > 0:
                        # NLOS transmission, function in DeLOSNLOS
                        Rec_power_busprv = NonLineofSight(Pt, busHeight, rxHeight, Gt, Gr, busprv_dist, walld_busprv, wallnum_busprv)
                    else:
                        # LOS transmission, function in DeLOSNLOS
                        Rec_power_busprv = LineofSight(Pt, busHeight, rxHeight, Gt, Gr, busprv_dist) 
                    if Rec_power_busprv > Rec_th:
                        busNeigh.append(veh)
                try:
                    # busNeigh didnt include the bus itself, but intNeigh includes.
                    HiTer = len(list(filter(lambda x: x not in busNeigh, intNeigh))) - 1
                except KeyError:
                    HiTer = 0
                Prob_RP = coll_prob(len(busNeigh), HiTer, T_pk)
                busls.append((curr_time, bus, busloc, intersec, TxRx_dist, Rec_power, walld, wallnum, 1-Prob_RP, len(busNeigh), HiTer))
            else:
                Prob_RP = 0
                busls.append((curr_time, bus, busloc, intersec, TxRx_dist, Rec_power, walld, wallnum, 1-Prob_RP))

traci.close()
labels = ['simtime', 'busID', 'buspos', 'intpos', 'dist', 'recpow', 'walldist', 'wallnum', 'Tot_prob', 'Neighbours', 'Hidden_Terminals']
df = pd.DataFrame.from_records(busls, columns=labels)
print("--- {}s seconds ---".format(time.time() - start_time))
print("--- Number of vehicles in each step ---", Num_veh)
# print("--- Number of buses in each step ---", Num_bus)

with open('update3012dataV2I.pkl', 'wb') as fw:
    pickle.dump(df, fw)
