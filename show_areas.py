# -*- coding: utf-8 -*-
"""
Created on Mon Jun 28 12:21:14 2021

For a given receiver deployment show the usable and coverage areas.
If 'bellhop=1' use the TL map from Bellhop. Else assume cylindrical spread
program adjusted for SD i.e.
grid size - according to the area's grid
TL - two options for the TL map
    a. TL of San-Diego for typical SVP
    b. TL of San-Diego for iso SVP
TL map was calculated using Bellhop package
Coverage files -
either set manually the receivers' positions by variable POP or read results of
'find_areas.py'. Here the results for 5 and 10 receivers for typical SVP and
for 10 receivers for iso SVP are provided.
TL map and Coverage results files are the
    'SD_tl_data_50x40_iso' and
    'cover_data_50x40.plk'
Respectively

@author: Talmon
"""

import geopandas as gp
import pandas as pd
from shapely.geometry import Polygon
import numpy as np
import matplotlib.pyplot as plt
import itertools
import utils

plt.ioff()  # swich off ploting for contour calculation

'''
make grid size same as calculated by BELLHOP with matlab
and load TL data too'''
# define our new play ground for real world analysis
x_grid = 51
y_grid = 41
mesh_size = x_grid  # some shortcut to claculate the visibility matrix
x = np.linspace(0, x_grid-1, x_grid)
y = np.linspace(0, y_grid-1, y_grid)
X, Y = np.meshgrid(x, y, indexing='ij')  # take care to matrix index
bellhop = 1  # bellhop = 1 for bellhop TL. otherwise use cylindrical TL
if bellhop == 1:
    TL_df = pd.read_pickle('SD_tl_data_50x40_iso.pkl')  # SD TL DF
etha = 0.95


''''
if POP to be read from file, run the two following lines.
if POP to be inserted manually, uncomment the lines bellow

'''
Result_df = pd.read_pickle('cover_data_50x40_iso_10R.pkl')
POP = Result_df['pop'].loc[np.argmax(Result_df['coverage'])]
# '''3R'''
# POP = [35, 31, 31, 25, 28, 30]  # PDAD
# POP = [29, 22, 21, 22, 25, 15]  # center
# # '''5R'''
# POP = [31, 26, 35, 19, 25, 37, 29, 23, 29, 25]  # PDAD
# POP = [30, 22, 25, 25, 20, 22, 22, 16, 28, 16]  # center
# # '''10R'''
# POP = [30, 21, 29, 24, 26, 25, 23, 24, 21, 22, 20, 19, 21, 16, 24, 15, 27,
#        16, 29, 18]  # center
# POP = [37, 17, 22, 13, 35, 17, 27, 38, 10, 26, 21, 29, 34, 25, 23, 19, 28,
#        25, 30, 38]  # PDAD
x_c = POP[::2]
y_c = POP[1::2]
N = len(x_c)
r = 5 * np.ones(N)  # used when cylindrical spreading is used
# standard deployment - equilateral shape
compare_to_standard = False
if compare_to_standard:
    rec_range_ratio = 0.5  # assumin max receiving range 1000m, grid size 100 m
    fact = 10 * rec_range_ratio
    x_std = np.zeros(N)
    y_std = np.zeros(N)

    def equi_penta(N):
        x_st = np.zeros(N)
        y_st = np.zeros(N)
        for i in range(N):
            x_st[i], y_st[i] = (fact * np.cos(2*np.pi*i/N+np.pi/(2*N)),
                                fact * np.sin(2*np.pi*i/N+np.pi/(2*N)))
        return x_st + x_grid/2, y_st + y_grid/2

    x_std, y_std = equi_penta(N)
    x_c = x_std.astype(int)
    y_c = y_std.astype(int)

# get all receiver combinatios
receivers_subset = ()
for i in range(3, N + 1):
    for subset in itertools.combinations(range(N), i):
        receivers_subset = receivers_subset + (subset, )

# get the usable_area area by finding intersections of receivers
receiver_df = gp.GeoDataFrame({'geometry': []})
for i in range(N):
    polys = utils.get_tl(x_c[i], y_c[i], r[i], bellhop, TL_df)
    receiver_df = receiver_df.append(gp.GeoDataFrame({'geometry': polys}),
                                     ignore_index=True)
p = Polygon([(0, 0), (1000, 0), (1000, 1000), (0, 1000)])
polys_aoi = gp.GeoDataFrame({'geometry': [p]})
usable_df = gp.GeoDataFrame({'geometry': []})
usable_df.insert(1, 'receiver set', '')
usable_df.insert(2, 'x coord', '')
usable_df.insert(3, 'y coord', '')
usable_df.insert(4, 'rounded x', '')
usable_df.insert(5, 'rounded y', '')
for k in range(len(receivers_subset)):
    usables = polys_aoi
    for i in receivers_subset[k]:
        usables = usables.intersection(receiver_df['geometry'].loc[i])
    usable_df = usable_df.append(gp.GeoDataFrame({'geometry': usables}),
                                 ignore_index=True)
    usable_df.at[k, 'receiver set'] = receivers_subset[k]
usable_df = usable_df.drop(usable_df[usable_df.is_empty].index.values,
                           axis=0).reset_index(drop=True)
usable_df = usable_df.loc[usable_df.geom_type == 'Polygon'
                          ].reset_index(drop=True)
num_of_polys = usable_df.shape[0]  # number of non-empty sub-usable areas
if num_of_polys == 0:
    print(' no usable area found')

GDOP_df = gp.GeoDataFrame({'geometry': []})
gdop_threshold = 5
for i in range(num_of_polys):
    x_xt, y_xt = utils.get_exterior(gp.GeoDataFrame(
        {'geometry': usable_df.loc[i]}), 'usable')
    usable_df.at[i, 'x coord'] = x_xt  # just for controling the rounded coord
    usable_df.at[i, 'y coord'] = y_xt
    usable_df.at[i, 'rounded x'] = np.array(x_xt).round(0)
    usable_df.at[i, 'rounded y'] = np.array(y_xt).round(0)
    x_vis = list(map(lambda k: x_c[k], usable_df['receiver set'].loc[i]))
    y_vis = list(map(lambda k: y_c[k], usable_df['receiver set'].loc[i]))
    H = utils.vis_mat(len(usable_df['receiver set'].loc[i]), x_vis, y_vis,
                      x_grid, y_grid, X, Y)
    x_xt_round = usable_df['rounded x'].loc[i]
    y_xt_round = usable_df['rounded y'].loc[i]
    # get coordinates in an ordered manner for GDOP colculation
    index_sets = [np.argwhere(k == y_xt_round) for k in np.unique(y_xt_round)]
    usable_area = np.zeros((len(index_sets), 3))
    for k in range(len(index_sets)):
        '''BELLHOP is not aware of the grid, thus some correction to stay
            inside the gride maybe required'''
        usable_area[k, 0] = y_xt_round[index_sets[k][0]]
        if y_xt_round[index_sets[k][0]] > y_grid - 1:  # stay inside the grid
            usable_area[k, 0] = y_grid - 1
        usable_area[k, 1] = min(x_xt_round[index_sets[k]])
        if min(x_xt_round[index_sets[k]]) > x_grid - 1:  # stay inside the grid
            usable_area[k, 1] = x_grid - 1
        usable_area[k, 2] = max(x_xt_round[index_sets[k]])
        if max(x_xt_round[index_sets[k]]) > x_grid - 1:  # stay inside the grid
            usable_area[k, 2] = x_grid - 1
        usable_area = usable_area.astype(int)
    GDOP = utils.gdop_map(H, usable_area, x_grid, y_grid)
    GDOP = np.nan_to_num(GDOP, nan=12, posinf=12, neginf=12)
    polys = utils.get_contour(X, Y, GDOP, gdop_threshold)
    if polys.is_valid[0]:
        GDOP_df = GDOP_df.append(gp.GeoDataFrame({'geometry': polys}),
                                 ignore_index=True)
# get the union of the areas
comp_usable_df = usable_df.drop(['receiver set', 'x coord', 'y coord',
                                 'rounded x', 'rounded y'], axis=1)
comp_usable_df = gp.GeoDataFrame({'geometry':
                                  gp.GeoSeries(comp_usable_df[
                                      'geometry'].unary_union)})
comp_usable_df = comp_usable_df.explode()
comp_usable_df = comp_usable_df.loc[comp_usable_df.geom_type == 'Polygon'
                                    ].reset_index(drop=True)
num_of_usables = comp_usable_df.shape[0]  # for plotting control
comp_GDOP_df = gp.GeoDataFrame({'geometry':
                                gp.GeoSeries(GDOP_df['geometry'].unary_union)})
comp_GDOP_df = comp_GDOP_df.explode()
comp_GDOP_df = comp_GDOP_df.loc[comp_GDOP_df.geom_type == 'Polygon'
                                ].reset_index(drop=True)
num_of_GDOP = comp_GDOP_df.shape[0]  # for plotting control
a_o_i = gp.GeoSeries(Polygon([(10, 10), (40, 10), (40, 40), (10, 40)]))
area_of_interest_df = gp.GeoDataFrame({'geometry': a_o_i})
if comp_usable_df.empty:
    comp_usable_df = utils.avoid_ter()
if comp_usable_df['geometry'].isna()[0]:
    comp_usable_df = utils.avoid_ter()
if comp_GDOP_df.empty:
    comp_GDOP_df = utils.avoid_ter()
if comp_GDOP_df['geometry'].isna()[0]:
    comp_GDOP_df = utils.avoid_ter()
penalty = utils.objective_penalty(comp_usable_df, comp_GDOP_df,
                                  area_of_interest_df, etha)
if not penalty:  # this happens when coverage area to aoi ratio > etha
    penalty = 0

# some ploting
plt.figure()
for i in range(N):
    x_xt, y_xt = utils.get_exterior(receiver_df['geometry'].loc[i], 'receiver')
    plt.plot(x_xt, y_xt)
x_xt, y_xt = utils.get_exterior(area_of_interest_df, 'aoi')  # exterior a_o_i
plt.plot(x_xt, y_xt)
for i in range(num_of_polys):
    plt.plot(usable_df['x coord'].loc[i], usable_df['y coord'].loc[i])
plt.plot(x_c, y_c, '.', label='Receiver positions')
idx = range(1, N+1)
for i, txt in enumerate(idx):
    plt.annotate(txt, (x_c[i], y_c[i]))
ax = plt.gca()
ax.set_aspect('equal')
plt.title('receiving areas intersection')
plt.show()

plt.figure()  # individual usable areas
for i in range(num_of_polys):
    # x_xt, y_xt = get_exterior(comp_usable_df['geometry'].loc[i], 'usable')
    plt.plot(usable_df['rounded x'].loc[i], usable_df['rounded y'].loc[i])
plt.title('usable areas')
plt.plot(x_c, y_c, '.', label='Receiver positions')
idx = range(1, N+1)
for i, txt in enumerate(idx):
    plt.annotate(txt, (x_c[i], y_c[i]))
plt.show()

plt.figure()  # union of all usable areas
for i in range(num_of_usables):
    x_xt, y_xt = utils.get_exterior(gp.GeoDataFrame(
        {'geometry': comp_usable_df.loc[i]}), 'usable')
    plt.plot(x_xt, y_xt)
plt.title('complete usable area contour')
plt.show()

plt.figure()  # usable and GDOP area
for i in range(num_of_usables):
    x_xt, y_xt = utils.get_exterior(gp.GeoDataFrame(
        {'geometry': comp_usable_df.loc[i]}), 'usable')  # exterior of usables
    plt.plot(x_xt, y_xt, 'g', label='Usable area {}'.format(i))
for i in range(num_of_GDOP):
    x_xt, y_xt = utils.get_exterior(gp.GeoDataFrame(
        {'geometry': comp_GDOP_df.loc[i]}), 'usable')  # exterior of GDOP areas
    plt.plot(x_xt, y_xt, 'b', label='GDOP = {}, area {}'
             .format(gdop_threshold, i))
x_xt, y_xt = utils.get_exterior(area_of_interest_df, 'aoi')  # exterior a_o_i
plt.plot(x_xt, y_xt)
plt.title('usable and GDOP areas contour')
plt.plot(x_c, y_c, '.', label='Receiver positions')
idx = range(1, N+1)
for i, txt in enumerate(idx):
    plt.annotate(txt, (x_c[i], y_c[i]))
# plt.plot(np.append(x_std, x_std[0]), np.append(y_std, y_std[0]))
# plt.axis('equal')
plt.legend()
plt.show()
# plt.savefig('standard method.pdf')

print('usable area', sum(comp_usable_df.area[:]))
print('usable to area of interest ration',
      sum(comp_usable_df.area[:]) / area_of_interest_df.area[0])
print('coverage area', sum(comp_GDOP_df.area[:]))
print('coverage to usability ratio',
      sum(comp_GDOP_df.area[:]) / sum(comp_usable_df.area[:]))
print('coverage to area of interest ratio',
      sum(comp_GDOP_df.area[:]) / area_of_interest_df.area[0])
print('penalty', penalty)
print('covrage inside a.o.i',  sum(gp.overlay(area_of_interest_df,
                                              comp_GDOP_df,
                                              how='intersection').area[:]))
print('receivers coordinates x', x_c)
print('receivers coordinates y', y_c)
