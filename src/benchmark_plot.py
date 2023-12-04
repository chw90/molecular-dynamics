#!/usr/bin/env python3

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from typing import List, Tuple
from scipy.stats import linregress

base = "argon"
extension = "log"


class LogFile:
    """Store metadata of a simulation log file"""

    path = None
    base = None
    extension = None
    data = pd.DataFrame()

    def __init__(self, path):
        self.path = path
        self.base = path.split("_")[0]
        self.extension = path.split(".")[-1]
        self.__read_data__()

    def __read_data__(self):
        self.data = pd.read_csv(self.path)
        self.data.columns = self.data.columns.str.strip()


def conforming(filename: str) -> bool:
    """Check if filename has the required base name and extension"""
    conforms = (
        filename.split("_")[0] == base
        and filename.split(".")[-1] == extension
        and filename.split("_")[-1].split(".")[0].isdigit()
    )
    return conforms


def get_particle_count_and_mean(lf: LogFile) -> Tuple[int, float]:
    particle_count = int(lf.path.split("_")[-1].split(".")[0])
    try:
        mean = lf.data["pair runtime"].mean()
    except KeyError:
        print("Log file " + lf.path + " has no 'pair runtime' column!")
        raise
    return (particle_count, mean)


def get_plot_data(dumps: List[LogFile]):
    ContainerVector = []
    ContainerCells = []
    for dump in dumps:
        data_tuple = get_particle_count_and_mean(dump)
        if dump.path.split("_")[1] == "cell":
            ContainerCells.append(data_tuple)
        else:
            ContainerVector.append(data_tuple)
    return (np.array(ContainerVector), np.array(ContainerCells))


def set_plot(dumps):
    data = get_plot_data(dumps)

    # sort for increasing particle count
    ContainerVector = np.sort(data[0], 0)
    ContainerCells = np.sort(data[1], 0)

    # plot data
    fig, ax = plt.subplots(1, 2)
    if ContainerVector.size > 0:
        xy = (ContainerVector[:, 0], ContainerVector[:, 1])
        slope = linregress(np.log(xy[0]), np.log(xy[1])).slope
        opt = {
            "color": "red",
            "label": rf"ContainerVector (slope $\approx$ {slope:4.3f})",
        }
        ax[0].plot(*xy, **opt)
        ax[1].loglog(*xy, **opt)
    if ContainerCells.size > 0:
        xy = (ContainerCells[:, 0], ContainerCells[:, 1])
        slope = linregress(np.log(xy[0]), np.log(xy[1])).slope
        opt = {
            "color": "green",
            "label": rf"ContainerCells (slope $\approx$ {slope:4.3f})",
        }
        ax[0].plot(*xy, **opt)
        ax[1].loglog(*xy, **opt)

    # add labels and legend
    ax[0].set_xlabel("particle count")
    ax[1].set_xlabel("logarithm of particle count")
    ax[0].set_ylabel("runtime")
    ax[1].set_ylabel("logarithm of runtime")
    ax[0].legend()
    ax[1].legend()

    # add grid for loglog plot
    ax[1].grid(which="both")
    return (fig, ax)


# get all dumps with the given base name and extension in current directory
dumps = [
    LogFile(filename) for filename in os.listdir(os.curdir) if conforming(filename)
]

fig, ax = set_plot(dumps)
plt.show()
