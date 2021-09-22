#!/usr/bin/env python

'''
Author:
    Shohei Kojima @ RIKEN

Requirement:
    Python 3.6 or later.

Description:
    This is an example script to plot histogram for similarity scores of LTR/ERVK.
    Example arbitrary threshold will be shown as a vertical dot line.

Usage:
    python %prog

License:
    See LICENSE file.
'''


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

threshold=200

# read .bed.gz file
f='file_name.bed.gz'
df=pd.read_csv(f, sep='\t', comment='#', header=None)
df.columns=['chr', 'start', 'end', 'name', 'score', 'strand']
df['TEclass']=[ i.split('.')[2] for i in df['name'] ]
df=df[ df['TEclass'] == 'LTR/ERVK' ]

# plot scores
ax=sns.histplot(x='score', data=df)
ax.axvline(x=threshold, linestyle='dashed')
ax.set_xlim(left=0, right=2000)
plt.savefig('hist.pdf')
