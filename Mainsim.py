"""This script provides the main simulation of SUMO and a proposed module. Note that for the purpose of comparison, 
the vehicle positions of the outputs are adjusted according to Veins."""
import os, sys, time, subprocess, pickle  
sys.path.append(os.path.join('c://','Program Files','sumo-0.30.0','tools'))
import traci
import math
from DeLOSNLOS import * 
from total_probnew import coll_prob
import time, random
from rtree import index
import pandas as pd
import itertools
"""Update version, includes rec_pow, positions, collision probability, and etc. 
   NOTE-1: the position of vehicles has been changed based on veins rule-- add margin, and y flipflop.
   NOTE-2: the rec_pow function has been modified, thus comment out the antenna gain part.
   NOTE-3: the dist_thre has been modified, based on two ray interference model.
   NOTE-4: move the dist_thre lines to up, in order to speed up a bit.
   NOTE-5: to speed up the model, add neighbours to each vehicle before calculate the rec.pow.
   NOTE-6: the x-y coordinates should be the origin version, but change them to fit veins positions for output.
   NOTE-7: updated building geometry file, the previous version missed a lot of building geometries.
   NOTE-8: to do sensitivity analysis, changes the vehicle density. 5 different simulation time, each time with 10 seconds.
   NOTE-9: fix the niehgbours and hidden terminals part, add limitation of receiving power level."""


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
# antenna height, in m
txHeight = 1.895
rxHeight = 1.895
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
sumoConfig = "path/to/sumocfg/file"
sumoProcess = subprocess.Popen([sumoBinary, '-c', sumoConfig, "--remote-port", str(PORT)], stdout= sys.stdout, stderr=sys.stderr)
# setup traci
traci.init(PORT)


# initialize an empty list
List = []
# initialize an empty list to contain number of vehicles in each step
Num_veh = []

for step in range(5010):
    # warm up sumo simulation
    traci.simulationStep()
    if step < 5000:
        continue
    if step == 5001:
    # record the time of simulation 
        start_time = time.time()
    # initialize an empty list to contain tuple of vehIDs that have enough receiving power level
    vehtup = []
    # initialize an empty list to contain tuple of vehIDs that have lower receiving power level
    failtup = []
    # initialize an empty dictionary contain receving power and vehIDs tuple
    recpower = {}
    default_xy = {}
    ID_Neighbours = {}
    # get all vehicle id list
    vehIDlist = traci.vehicle.getIDList()
    # sumo simulation time
    curr_time = traci.simulation.getCurrentTime()
    for vehid in vehIDlist:
        # get vehicle x-y position, and write in a dictionary, add Veins' werid margin value
        x, y = traci.vehicle.getPosition(vehid)
        # check if the vehicle's position is in the region of interest area (use true values without margin)
        if meetsROIareaxml((x, y), ROI):
            default_xy.update({vehid:(x, y)})
    xypositions = list(default_xy.values()) 
    # number of vehicles in the ROI of each step
    Num_veh.append(len(xypositions))

    # check all possible communication pairs and write them into a list of tuples, consider both transmission range and receivng power level
    for (index, index_2) in list(itertools.combinations(range(len(xypositions)), 2)):
        sendpos = xypositions[index]
        recpos = xypositions[index_2]
        # calculate the distance between transmitting and receving nodes, in meters, function in DeLOSNLOS
        TxRx_dist = caldistance(sendpos, recpos)
        # ignore comm_pair with distance>977, in this way also speed up model and the results would be correct.
        if TxRx_dist > dist_th:
            continue
        # NOTE: might speed up the model.-- after checking in my computer, only speed up about 1 second than commented out lines.
        vehID = [k for k,v in default_xy.items() if v == sendpos][0]
        vehID_2 = [k for k,v in default_xy.items() if v == recpos][0]
        # Get the total length of intersection and the number of crossing walls of the building, x-y coordinates
        walld, wallnum = LOSNLOSxml(js, idx, sendpos, recpos)
        # Based on the transmission way, calculate the receiving power level
        if wallnum > 0:
            # NLOS transmission, function in DeLOSNLOS
            Rec_power = NonLineofSight(Pt,txHeight, rxHeight, Gt, Gr, TxRx_dist, walld, wallnum)
        else:
            # LOS transmission, function in DeLOSNLOS
            Rec_power = LineofSight(Pt, txHeight, rxHeight, Gt, Gr, TxRx_dist) 
        if Rec_power > Rec_th:
            # vehtup includes the communicating pair that has a higher receiving power lvl than the threshold
            vehtup.append((vehID, vehID_2))
            recpower.setdefault((vehID, vehID_2), (Rec_power, TxRx_dist, walld, wallnum))
            recpower.setdefault((vehID_2, vehID), (Rec_power, TxRx_dist, walld, wallnum))
        else:
            failtup.append((vehID, vehID_2))
            recpower.setdefault((vehID, vehID_2), (Rec_power, TxRx_dist, walld, wallnum))
            recpower.setdefault((vehID_2, vehID), (Rec_power, TxRx_dist, walld, wallnum))
    # get all communicating vehicle's (rec.power > rec.threshold) unique ID in list survID
    slist = [item for sublist in vehtup for item in sublist]
    survID = list(set(slist))
    # create an empty dictionary for all the vehicles' neighbours
    Neighbours = {}
    # get all neighbours of survID, the neighbours that are able to hear from the survID
    for ids in survID:
        # neighpair: list of tuples that include the ids
        neighpair = [item for item in vehtup if ids in item]
        # unique neighbours of ids
        neigh = list(set([item for sublist in neighpair for item in sublist]))
        Neighbours.setdefault(ids, neigh)
    for tran, recvls in Neighbours.items():
        sendpos = default_xy[tran]
        sendx = round((sendpos[0] + 25), 2)
        sendy = round((11455 - sendpos[1] + 25), 2)
        Sendpos = (sendx, sendy)
        tranneigh = len(recvls)
        for recID in recvls:
            if recID == tran:
                continue
            recpos = default_xy[recID]
            recx = round((recpos[0] + 25), 2)
            recy = round((11455 - recpos[1] + 25), 2)
            Recpos = (recx, recy)
            recvneigh = Neighbours[recID]
            try:
                HiTer = len(list(filter(lambda x: x not in recvls, recvneigh)))
            except KeyError:
                HiTer = 0
            power = recpower[(tran, recID)][0]
            dist = recpower[(tran, recID)][1]
            wd = recpower[(tran, recID)][2]
            wn = recpower[(tran, recID)][3]

            Prob_RP = coll_prob(tranneigh, HiTer, T_pk)
            # generate a random number between 0 and 1
            p = random.uniform(0, 1)
            # based on random number and prob. decide if this trial is successful or not
            randres = bool (p < Prob_RP)
            List.append((curr_time, tran, recID, power, Sendpos, Recpos, dist, wd, wn, Prob_RP, tranneigh, HiTer, randres))


traci.close()
labels = ['simtime', 'sendID', 'recID', 'recpow', 'sendpos', 'recpos', 'dist', 'walldist', 'wallnum', 'Prob_RP', 'Neighbours', 'Hidden_Terminals', 'randresult']
df = pd.DataFrame.from_records(List, columns=labels)
print("--- {}s seconds ---".format(time.time() - start_time))
print("--- Number of vehicles in each step ---", Num_veh)

with open('output pickle', 'wb') as fw:
    pickle.dump(df, fw)
