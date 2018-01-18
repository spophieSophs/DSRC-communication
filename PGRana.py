import pickle
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
"""This script focus on the impacts of the packet generation rate on the total collision probability."""

if __name__ == '__main__':
    with open('update3012dataV2I.pkl', 'rb') as f10, open("5Hz1001data.pkl", 'rb') as f5, open('2Hz1001data.pkl', 'rb') as f2, open('1Hz1001data.pkl', 'rb') as f1:
        df10Hz = pickle.load(f10)
        df5Hz = pickle.load(f5)
        df2Hz = pickle.load(f2)
        df1Hz = pickle.load(f1)
    df10Hz = df10Hz[df10Hz.Tot_prob != 1]
    df5Hz = df5Hz[df5Hz.Tot_prob != 1]
    df2Hz = df2Hz[df2Hz.Tot_prob != 1]
    df1Hz = df1Hz[df1Hz.Tot_prob != 1]
    
 
    Hz = ['10Hz', '5Hz', '2Hz', '1Hz']
    intA10Hz = df10Hz[df10Hz.intpos == (6940, 6902)]
    intA5Hz= df5Hz[df5Hz.intpos == (6940, 6902)]
    intA2Hz = df2Hz[df2Hz.intpos == (6940, 6902)]
    intA1Hz = df1Hz[df1Hz.intpos == (6940, 6902)]
    intAHz = pd.concat([intA10Hz['Tot_prob'], intA5Hz['Tot_prob'], intA2Hz['Tot_prob'], intA1Hz['Tot_prob']], axis=1, keys =Hz)
    intB10Hz = df10Hz[df10Hz.intpos == (7074, 7257)]
    intB5Hz = df5Hz[df5Hz.intpos == (7074, 7257)]
    intB2Hz = df2Hz[df2Hz.intpos == (7074, 7257)]
    intB1Hz = df1Hz[df1Hz.intpos == (7074, 7257)]
    intBHz = pd.concat([intB10Hz['Tot_prob'], intB5Hz['Tot_prob'], intB2Hz['Tot_prob'], intB1Hz['Tot_prob']], axis=1, keys =Hz)
    intC10Hz = df10Hz[df10Hz.intpos == (6971, 7602)]
    intC5Hz = df5Hz[df5Hz.intpos == (6971, 7602)]
    intC2Hz = df2Hz[df2Hz.intpos == (6971, 7602)]
    intC1Hz = df1Hz[df1Hz.intpos == (6971, 7602)]
    intCHz = pd.concat([intC10Hz['Tot_prob'], intC5Hz['Tot_prob'], intC2Hz['Tot_prob'], intC1Hz['Tot_prob']], axis=1, keys =Hz)
    
    name = ['A', 'B', 'C']
    fig = plt.figure(figsize=(7, 5))
    ax0 = fig.add_subplot(2, 2, 1) 
    ax0 = intAHz.boxplot()
    ax0.set_yticks(np.arange(0, 0.24, 0.04))
    ax0.set_ylim((0, 0.22))
    # ax0.set_title('Intersection A')
    ax0.set_xlabel('Frequency')
    ax0.set_ylabel('Total Collision Probability')
    ax1 = fig.add_subplot(2, 2, 2)
    ax1 = intBHz.boxplot()
    ax1.set_yticks(np.arange(0, 0.24, 0.04))
    ax1.set_ylim((0, 0.22))
    ax1.set_xlabel('Frequency')
    ax1.set_ylabel('Total Collision Probability')
    # ax1.set_title('Intersection B')
    ax2 = fig.add_subplot(2, 2, 3)
    ax2 = intCHz.boxplot()
    ax2.set_yticks(np.arange(0, 0.24, 0.04))
    ax2.set_ylim((0, 0.22))
    ax2.set_xlabel('Frequency')
    ax2.set_ylabel('Total Collision Probability')
    # ax2.set_title('Intersection C')
 
    
    # # ax.set_xlim((0, 200))
    # 
    # axls[1].set_title('Distance VS. Total Collision Probability', fontsize=10)
    plt.grid(True)
    # plt.xlabel('Packet Generation Frequency')
    # plt.ylabel('Total Collision Probability')
    plt.legend()
    fig.tight_layout()
    fig.savefig('1701PGRboxplot.png', dpi=300)
    plt.show()
 