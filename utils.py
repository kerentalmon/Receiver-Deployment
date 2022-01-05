# -*- coding: utf-8 -*-
"""
Created on Wed Jan  5 11:38:24 2022
helper function to find GDOP map, get the TL and the areas exteriors etc.
@author: Talmon
"""
import numpy as np
import geopandas as gp
from shapely.geometry import Polygon
import matplotlib.pyplot as plt

def get_tl(a, b, m, bellhop, TL_df):
    if bellhop == 1:
        # polys_x = TL_df['TL contour_x'].loc[a, b]
        # polys_y = TL_df['TL contour_y'].loc[a, b]
        # P = map(lambda k: (polys_x[0][k], polys_y[0][k]),
        #         range(np.size(polys_x)))
        ''' an alternative to the above - same result'''
        P = zip(TL_df['TL contour_x'].loc[a, b].squeeze(),
                TL_df['TL contour_y'].loc[a, b].squeeze())
        polys = gp.GeoSeries([Polygon(P)])
    else:
        F = circ(a, b, m)
        polys = get_contour(F, 0)
    return polys


def circ(a, b, r, X, Y):
    return (X - a) ** 2 + (Y - b) ** 2 - r ** 2


# get the exterior line of a polygon
def get_exterior(df, string):
    try:
        x = df['geometry'].iloc[0].exterior.coords.xy[0]
        y = df['geometry'].iloc[0].exterior.coords.xy[1]
    except Exception:
        print('no {} area' .format(string))
        x = np.array([0, 1, 1, 0])
        y = np.array([0, 0, 1, 1])
    return x, y


def get_contour(X, Y, C, threshold):
    try:
        p = plt.contour(X, Y, C, [threshold]).allsegs[0][0]
        plt.close()  # otherwise "contour" plots a plot
    except (UserWarning, IndexError, Exception):
        p = [(0, 0), (1, 0), (1, 1), (0, 1)]  # to avoid program termination
    if len(p) < 3:  # check that there is a GDOP area/ receiving area
        if threshold == 0:
            print('no receiver defined')
        else:
            print('no GDOP found')
        p = [(0, 0), (1, 0), (1, 1), (0, 1)]  # to avoid program termination
    polys = gp.GeoSeries([Polygon(p)])
    return polys


def avoid_ter():
    p = Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])
    polys = gp.GeoSeries([Polygon(p)])
    return gp.GeoDataFrame({'geometry': polys})


def vis_mat(n, x_r, y_r, x_grid, y_grid, X, Y):
    range_to_recevir = np.zeros([n, x_grid, y_grid])
    h_x = np.zeros([n, x_grid, y_grid])
    h_y = np.zeros([n, x_grid, y_grid])
    h = np.zeros([n, 2, x_grid, y_grid])
    for i in range(n):
        range_to_recevir[i, :, :] = np.sqrt((X - x_r[i])**2 + (Y-y_r[i])**2)
        inx = np.where(range_to_recevir[i, :, :] < 10e-4)
        range_to_recevir[i, inx] = 10e-4  # avoids division by zero
        h_x[i, :, :] = (X - x_r[i]) / range_to_recevir[i, :, :]
        h_y[i, :, :] = (Y - y_r[i]) / range_to_recevir[i, :, :]
        h[i, 0, :, :] = h_x[i, :, :]  # can omit h_x, h_y, here for clarity
        h[i, 1, :, :] = h_y[i, :, :]
    return h


def gdop_map(H, usable_area, x_grid, y_grid):
    GD = np.zeros([x_grid, y_grid]) + 17  # 17 is just for distinction
    for count in range(len(usable_area)):
        j = usable_area[count, 0]
        stp = 1
        if usable_area[count, 1] > usable_area[count, 2]:
            stp = -1
        for k in range(usable_area[count, 1], usable_area[count, 2]+stp, stp):
            try:
                D = np.linalg.inv(np.matmul(H[:, :, k, j].T, H[:, :, k, j]))
                GD[k, j] = np.sqrt(D[0, 0] ** 2 + D[1, 1] ** 2)
            except np.linalg.LinAlgError:
                GD[k, j] = 17  # 17 is just for distinction
    return GD


def objective_penalty(comp_usable_df, comp_GDOP_df, area_of_interest_df, etha):
    global penalty_0, penalty_00
    ''' penalty to force covrage area to be inside the area of interest'''
    penalty = 0
    penalty_0 = 0
    penalty_1 = 0
    coverage_area = sum(comp_GDOP_df.area[:])
    if coverage_area / area_of_interest_df.area[0] < etha:
        if not gp.overlay(area_of_interest_df, comp_usable_df,
                          how='intersection', keep_geom_type=False).empty:
            penalty_0 = sum(comp_usable_df.area[:] -
                            gp.overlay(area_of_interest_df, comp_usable_df,
                            how='intersection').area[:])
        else:
            penalty_0 = sum(comp_usable_df.area[:])
        if not gp.overlay(comp_GDOP_df, area_of_interest_df,
                          how='intersection', keep_geom_type=False).empty:
            penalty_1 = coverage_area -\
                sum(gp.overlay(area_of_interest_df, comp_GDOP_df,
                               how='intersection').area[:])
        else:
            penalty_1 = coverage_area
        penalty = penalty_0 + penalty_1
        return penalty
