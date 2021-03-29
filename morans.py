
import esda

import geopandas as gpd

import libpysal as lps


def process_adj(file):
    gol_ad = gpd.read_file(file)
    gol_ad = gpd.GeoDataFrame(gol_ad, geometry='geometry', crs={'init': 'epsg:2180'})
    listo = [72, 105, 109, 123, 127, 246, 374, 601]
    gol_ad = gol_ad.loc[~gol_ad.index.isin(listo)]
    we_ad = lps.weights.DistanceBand.from_dataframe(gol_ad, threshold=100, binary=False, alpha=-10)
    lm_ad = esda.Moran_Local(gol_ad['adjacency'], we_ad, transformation="r", permutations=99)
    lag_adj = lps.weights.lag_spatial(we_ad, gol_ad['adjacency'])

    sig = lm_ad.p_sim < 0.05

    coldspot = sig * lm_ad.q == 3


    spotsa = ['n.sig.', 'cold spot']
    labelsa = [spotsa[i] for i in coldspot * 1]
    gol_ad = gol_ad
    #gol_ad = gol_ad.to_crs(epsg=3857)
    # hmap = colors.ListedColormap(['blue', 'lightgrey'])
    # f, ax = plt.subplots(1, figsize=(9, 9))
    gol_ad = gol_ad.assign(cl=labelsa)  # .plot(column='cl', categorical=True, \
    # k=2, cmap=hmap, linewidth=0.1, ax=ax, \
    # edgecolor='white', legend=True, alpha=0.5)
    # ctx.add_basemap(ax, zoom=6)
    # ax.set_axis_off()
    # ctx.add_basemap(ax, zoom=15)
    gol_ad = gol_ad.loc[gol_ad['cl'] == 'cold spot']
    return gol_ad


def clust_age(file):
    gdf = gpd.read_file(file)
    gdf = gpd.GeoDataFrame(gdf, geometry='geometry', crs={'init': 'epsg:2180'})
    sh = gdf[['x', 'y', 'year_built', 'geometry']]
    # sh = gpd.GeoDataFrame(sh, geometry=gpd.points_from_xy(sh['y'], sh['x'],crs="EPSG:2180"))
    # sh = sh.set_crs(epsg=2180)

    sh = sh.loc[sh['year_built'] > 1850]
    sh['year_built'] = sh['year_built'].astype(int)
    #listo = [1120, 1442, 2429, 2442, 2506, 2558]
    listo = [5,361,750]
    sh = sh.loc[~sh.index.isin(listo)]
    we = lps.weights.DistanceBand.from_dataframe(sh, threshold=150, binary=False, alpha=-1.5)
    lm = esda.Moran_Local(sh['year_built'], we, transformation="r", permutations=99)
    lag_year = lps.weights.lag_spatial(we, sh['year_built'])

    sig = lm.p_sim < 0.05

    coldspot = sig * lm.q == 3

    spots = ['n.sig.', 'cold spot']
    labels = [spots[i] for i in coldspot * 1]
    # sh = sh
    ##sh = sh.to_crs(epsg=3857)
    # hmap = colors.ListedColormap(['blue', 'lightgrey'])
    # f, ax = plt.subplots(1, figsize=(9, 9))
    nf = sh.assign(cl=labels)  # .plot(column='cl', categorical=True, \
    # k=2, cmap=hmap, linewidth=0.1, ax=ax, \
    # edgecolor='white', legend=True, alpha=0.5)
    # ctx.add_basemap(ax, zoom=6)
    # ax.set_axis_off()
    # ctx.add_basemap(ax, zoom=15)
    # plt.show()
    nf = nf.loc[nf['cl'] == 'cold spot']

    return nf


def buffer(gdf):
    new_gdf = gdf.copy()
    new_gdf['geometry'] = new_gdf.geometry.buffer(50)
    mp = new_gdf.unary_union
    return mp


def calcarea(mg):
    df = gpd.GeoDataFrame(mg, crs="EPSG:2180", geometry=0)
    df['trace'] = 'd'
    df1 = df.dissolve(by='trace')
    return df1


def calcarea2(gdf1, gdf2):
    centers = gpd.clip(gdf1, gdf2)
    centers = centers.explode()
    centers = gpd.GeoDataFrame(centers, crs="EPSG:2180", geometry=0)

    centers['AREA'] = centers.geometry.area
    centers = centers.loc[centers['AREA'] > 20000]
    return centers.AREA.sum().round(2)