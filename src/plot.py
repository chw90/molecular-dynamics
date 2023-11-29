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

    fig, axes = plt.subplots(nrows=2, ncols=3)
    df.plot(ax=axes[0, 0], x='time', y='total energy', ylabel='total energy', legend=False)
    df.plot(ax=axes[0, 1], x='time', y='kinetic energy', ylabel='kinetic energy', legend=False)
    df.plot(ax=axes[0, 2], x='time', y='potential energy', ylabel='potential energy', legend=False)
    df.plot(ax=axes[1, 0], x='time', y='temperature', ylabel='temperature', legend=False)
    df.plot(ax=axes[1, 1], x='time', y='pressure', ylabel='pressure', legend=False)
    df.plot(ax=axes[1, 2], x='time', y='volume', ylabel='volume', legend=False)

    if args.avg:
        if not args.window:
            mean = df.mean()

            color = plt.rcParams['axes.prop_cycle'].by_key()['color'][1]
            axes[0, 0].axhline(y=mean["total energy"], color=color)
            axes[0, 1].axhline(y=mean["kinetic energy"], color=color)
            axes[0, 2].axhline(y=mean["potential energy"], color=color)
            axes[1, 0].axhline(y=mean["temperature"], color=color)
            axes[1, 1].axhline(y=mean["pressure"], color=color)
            axes[1, 2].axhline(y=mean["volume"], color=color)
        else:
            mean = df.rolling(
                args.window,
                min_periods=(args.window // 2),
                on='step',
                center=True
            ).mean()

            mean.plot(ax=axes[0, 0], x='time', y='total energy', legend=False)
            mean.plot(ax=axes[0, 1], x='time', y='kinetic energy', legend=False)
            mean.plot(ax=axes[0, 2], x='time', y='potential energy', legend=False)
            mean.plot(ax=axes[1, 0], x='time', y='temperature', legend=False)
            mean.plot(ax=axes[1, 1], x='time', y='pressure', legend=False)
            mean.plot(ax=axes[1, 2], x='time', y='volume', legend=False)

    hsep = 0.05
    vsep = 0.1
    plt.subplots_adjust(left=hsep, right=1 - hsep, top=1 - vsep, bottom=vsep)
    plt.show()
