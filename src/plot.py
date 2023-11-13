#!/usr/bin/env python3

import argparse
import pandas as pd
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()

parser.add_argument(
    'filepath',
    help='path to the statistics output file'
)
parser.add_argument(
    '--avg',
    action='store_true',
    help='plot mean or if --window is set rolling mean'
)
parser.add_argument(
    '--window',
    type=int,
    help='window size for rolling mean'
)
args = parser.parse_args()

with open(args.filepath) as f:
    df = pd.read_csv(f)
    df.columns = df.columns.str.strip()

    fig, axes = plt.subplots(nrows=1, ncols=3)
    df.plot(ax=axes[0], x='time', y='kinetic energy', ylabel='Kinetic Energy', legend=False)
    df.plot(ax=axes[1], x='time', y='temperature', ylabel='Temperature', legend=False)
    df.plot(ax=axes[2], x='time', y='pressure', ylabel='Pressure', legend=False)

    if args.avg:
        if not args.window:
            mean = df.mean()

            color = plt.rcParams['axes.prop_cycle'].by_key()['color'][1]
            axes[0].axhline(y=mean["kinetic energy"], color=color)
            axes[1].axhline(y=mean["temperature"], color=color)
            axes[2].axhline(y=mean["pressure"], color=color)
        else:
            mean = df.rolling(
                args.window,
                min_periods=(args.window // 2),
                on='step',
                center=True
            ).mean()

            mean.plot(ax=axes[0], x='time', y='kinetic energy', legend=False)
            mean.plot(ax=axes[1], x='time', y='temperature', legend=False)
            mean.plot(ax=axes[2], x='time', y='pressure', legend=False)

    sep = 0.05
    plt.subplots_adjust(left=sep, right=1 - sep, top=1 - sep, bottom=sep)
    plt.show()
