#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 27 January 2025

@author: liciaray
"""

import pandas as pd

df = pd.read_csv('../input/2004_2012_cassini_position_in_one_hour_bins.csv')


df_short = df[['TIME','orb_x','orb_y','orb_z']]

df_short['TIME'] = pd.to_datetime(df_short['TIME'], format='%d/%m/%Y %H:%M')
#df_short['TIME'] = df_short['TIME'].dt.strftime('%Y-%m-%dT%H:%M')

#df_short['TIME'] = pd.to_datetime(df_short['TIME'], format='%Y-%m-%dT%H:%M:%S')
df_short['TIME'] = df_short['TIME'].dt.strftime('%Y-%m-%dT%H:%M')

numeric_columns = df_short.select_dtypes(include=['number']).columns
print(numeric_columns)

#df_short['orb_x'] = df_short['orb_x'].apply(lambda x: f"{x:.8f}")
#df_short['orb_y'] = df_short['orb_y'].apply(lambda x: f"{x:.8f}")
#df_short['orb_z'] = df_short['orb_z'].apply(lambda x: f"{x:.8f}")

df_short['orb_x'] = df_short['orb_x'].apply(lambda x: f"{x: >13.8f}")
df_short['orb_y'] = df_short['orb_y'].apply(lambda x: f"{x: >13.8f}")
df_short['orb_z'] = df_short['orb_z'].apply(lambda x: f"{x: >13.8f}")

df_short.to_csv('../input/list_for_tracing.csv', sep=',', index=False, header=True)





