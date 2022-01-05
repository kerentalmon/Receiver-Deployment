# -*- coding: utf-8 -*-
"""
Created on Mon Jun 28 12:16:38 2021
Optimal receiver deployment search for N receivers using evolutionary
algorithms base on DEAP and SCOOP packages
Adjusted to San-Diego thus, uses 50X40 grid with BELLHOP TL calculation
use the 'bellhop' flag to decide if TL map is used or cylinder spread in use
Take care to adjust:
    - TL file name
    - grid size
    - area of interest size
    - bound for evolutionary algo
    - number of receivers N
    - uncomment the future.map line of wishing to run on parallel cores
    - input ga pramenters e.g. number of generations
@author: Talmon
"""
import array
import random
import pandas as pd
from deap import algorithms
from deap import base
from deap import creator
from deap import tools
from scoop import futures  # used for server runs
import geopandas as gp
from shapely.geometry import Polygon
import numpy as np
import matplotlib.pyplot as plt
import time
import itertools
import warnings
import utils
warnings.filterwarnings('ignore')
try:
    del creator.FitnessMax
    del creator.Individual
except AttributeError:
    pass


def checkBounds(x_max, y_max):
    def decorator(func):
        def wrapper(*args, **kargs):
            offspring = func(*args, **kargs)
            for child in offspring:
                for i in range(2 * N):
                    if i < N:
                        if child[i] >= x_max:
                            child[i] = x_max - 1
                        elif child[i] < 0:
                            child[i] = 0
                    else:
                        if child[i] >= y_max:
                            child[i] = y_max - 1
                        elif child[i] < 0:
                            child[i] = 0
            return offspring
        return wrapper
    return decorator


def find_coverage(POP):
    to_append = list(ind for ind in POP)
    global df_1
    df_1.at[len(df_1), 'pop'] = to_append
    x_c = POP[::2]
    y_c = POP[1::2]
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
    # create dummy df for a start - the complete playground
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
                               axis=0)
    usable_df = usable_df.loc[usable_df.geom_type == 'Polygon'
                              ].reset_index(drop=True)
    num_of_polys = usable_df.shape[0]  # number of non-empty sub-usable areas
    # if num_of_polys == 0:
    #     print(' no usable area found')
    GDOP_df = gp.GeoDataFrame({'geometry': []})
    gdop_threshold = 5
    for i in range(num_of_polys):
        x_xt, y_xt = utils.get_exterior(gp.GeoDataFrame(
            {'geometry': usable_df.loc[i]}), 'usable')
        usable_df.at[i, 'x coord'] = x_xt  # for controling the rounded coord
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
        index_sets = [np.argwhere(k == y_xt_round) for k in
                      np.unique(y_xt_round)]
        usable_area = np.zeros((len(index_sets), 3))
        for k in range(len(index_sets)):
            '''BELLHOP is not aware of the grid, thus some correction to stay
            inside the gride maybe required'''
            # usable_area: from min x_xt to max x_xt at hight of y_xt - rounded
            usable_area[k, 0] = y_xt_round[index_sets[k][0]]
            if y_xt_round[index_sets[k][0]] > y_grid - 1:
                usable_area[k, 0] = y_grid - 1  # stay inside the grid
            usable_area[k, 1] = min(x_xt_round[index_sets[k]])
            if min(x_xt_round[index_sets[k]]) > x_grid - 1:
                usable_area[k, 1] = x_grid - 1  # stay inside the grid
            usable_area[k, 2] = max(x_xt_round[index_sets[k]])
            if max(x_xt_round[index_sets[k]]) > x_grid - 1:
                usable_area[k, 2] = x_grid - 1  # stay inside the grid
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
    # num_of_usables = comp_usable_df.shape[0]  # for plotting control
    comp_GDOP_df = gp.GeoDataFrame(
        {'geometry': gp.GeoSeries(GDOP_df['geometry'].unary_union)})
    comp_GDOP_df = comp_GDOP_df.explode()
    comp_GDOP_df = comp_GDOP_df.loc[comp_GDOP_df.geom_type == 'Polygon'
                                    ].reset_index(drop=True)
    # num_of_GDOP = comp_GDOP_df.shape[0]  # for plotting control
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
    coverage_area = sum(comp_GDOP_df.area[:]) - penalty
    df_1.at[len(df_1)-1, 'coverage'] = coverage_area
    return coverage_area,


plt.ioff()  # swich off ploting

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
    TL_df = pd.read_pickle('SD_tl_data_50x40_iso.pkl')  # reading the SY TL DF


N = 5  # set number of receivers
etha = 0.95
r = 5 * np.ones(N)
df_1 = pd.DataFrame(data=[], columns=list(range(2*N)))
df_1 = pd.DataFrame(data=[], columns=['pop', 'coverage'])
df_1['pop'] = df_1['pop'].astype('object')

a_o_i = gp.GeoSeries(Polygon([(10, 10), (40, 10), (40, 40), (10, 40)]))
area_of_interest_df = gp.GeoDataFrame({'geometry': a_o_i})

creator.create("FitnessMax", base.Fitness, weights=(1.0,))
creator.create("Individual", array.array, typecode='i',
               fitness=creator.FitnessMax)
toolbox = base.Toolbox()
toolbox.register("attr_int_x", random.randint, 5, 45)  # Attribute generator
toolbox.register("attr_int_y", random.randint, 5, 35)  # Attribute generator

# Structure initializers
toolbox.register("individual", tools.initCycle, creator.Individual,
                 (toolbox.attr_int_x, toolbox.attr_int_y), n=N)
toolbox.register("population", tools.initRepeat, list, toolbox.individual)
toolbox.register("evaluate", find_coverage)
toolbox.register("mate", tools.cxTwoPoint)
toolbox.register("mutate", tools.mutUniformInt, low=5, up=40, indpb=0.4)
toolbox.register("select", tools.selTournament, tournsize=3)
toolbox.decorate("individual", checkBounds(x_grid, y_grid))
#  make sure next generations are inside the area of interest
toolbox.decorate("mate", checkBounds(x_grid, y_grid))
toolbox.decorate("mutate", checkBounds(x_grid, y_grid))
population = toolbox.population(n=N)
# toolbox.register("map", futures.map)  # uncomment for SCOOP paralle computing


def main():
    start = time.time()
    random.seed(64)
    pop = toolbox.population(n=10)
    hof = tools.HallOfFame(1)
    stats = tools.Statistics(lambda ind: ind.fitness.values)
    stats.register("avg", np.mean)
    stats.register("std", np.std)
    stats.register("min", np.min)
    stats.register("max", np.max)
    pop, log = algorithms.eaSimple(pop, toolbox, cxpb=0.7, mutpb=0.7, ngen=4,
                                   stats=stats, halloffame=hof, verbose=True)
    stop = time.time()
    print('run time', stop - start)
    return pop, log, hof


if __name__ == "__main__":
    pop, log, hof = main()
    print(max(df_1['coverage']))
    print(np.argmax(df_1['coverage']))
    print(df_1.loc[np.argmax(df_1['coverage'])])
    print(hof[-1])
    # df_1.to_pickle('cover_data_50x40_iso.pkl')  # save the data to HD
# uncomment future.map to run with scoop and run the below for command line
# python -m scoop coverage_area_2.py # for the server
