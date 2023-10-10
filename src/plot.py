#!/usr/bin/env python3

import argparse
import pandas as pd
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

with open(args.filename) as f:
    df = pd.read_csv(f)
    df.columns = df.columns.str.strip()

    fig, axes = plt.subplots(nrows=1, ncols=3)
    df.plot(ax=axes[0], x='time', y='kinetic energy', ylabel='Kinetic Energy')
    df.plot(ax=axes[1], x='time', y='temperature', ylabel='Temperature')
    df.plot(ax=axes[2], x='time', y='pressure', ylabel='Pressure')

    plt.show()
