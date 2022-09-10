"""Helper functions for mapper analysis of plant RNAseq data.

This file contains helper functions used in the mapper analysis of the
plant RNAseq data.

Functions
---------
    - loaddata() : Loads and preprocesses factors and RNAseq data.
    - colorscale_from_matplotlib_cmap(): Generates a keppler mapper colorscale
        from matplotlib cmap.
    - get_color_dict(): Returns a dictionary mapping factor levels to colors    
        from the tab20 colormap.
    - create_pie_chart_json(): Creates json files required to visualize mapper 
        graphs with pie chart nodes.

"""
import os
import numpy as np
import pandas as pd
import pickle
import json
import matplotlib.pyplot as plt
from matplotlib.colors import to_hex, to_rgb
from sklearn.preprocessing import scale
from scipy import stats
from itertools import combinations


def loaddata(factor_file: str, rna_file: str, factor_list: list,
             logscale: bool = False, zscore: bool = True, group: bool = False) -> (object, list):
    """Load and Preprocess Factor and RNAseq Data.

    This function is used to load preprocess the data provided in two input
    files. First file contains factors data. Second file contains RNAseq data.
    The data points in two files are matched by a unique 'sra'.

    The function performs an inner join of two files on the 'sra'.
    This extracts records for which factors and RNAseq are both available.

    Then, preprocessing steps (scaling, normalization, grouping) are performed.
    A single dataframe, containing both factors and RNAseq is returned along
    with a list of orthogroup names.

    Parameters
    ----------
    factor_file : str
        full path to the csv file containing factor data
    rna_file : str
        full path to the csv file containint RNAseq data
    factor_list : list
        List of factor names to be kept (column names in factor_file)
    logscale : bool, optional
        Scale RNAseq data on a log2 scale (y = np.log2(x+1)).
        The default is True.
    zscore : bool, optional
        Standardize RNAseq data.
        Each 'sra' is normalized to 0 mean 1 standard deviation.
        The default is True.
    group : bool, optional
        Group records by sample_id.
        There's a one-to-many mapping between sample_id and sra.
        If true, factors data is grouped by sample id (keeping first record)
        RNAseq data is averaged across sra.
        The default is True.

    Returns
    -------
    df : pandas dataframe
        Contains factors and orthogroups data after preprocessing.
    ortho_names : list
        Contains list of names of orthogroups.
    """
    # Read the two files into pandas dataframes
    try:
        factor_input = pd.read_csv(factor_file)
    except FileNotFoundError:
        print('Could not find {}'.format(factor_file))
        return None, None
    else:
        fnr, fnc = factor_input.shape
        print('Factor input data shape: ({}, {})'.format(fnr, fnc))

    try:
        rna_input = pd.read_csv(rna_file)
    except FileNotFoundError:
        print('Could not find {}'.format(rna_file))
        return None, None
    else:
        rnr, rnc = rna_input.shape
        print('RNA input data shape: ({}, {})'.format(rnr, rnc))

    # Transpose the RNAseq data so that rows are indexed by 'sra'
    rna_input = rna_input.set_index('Orthogroup').T
    rna_input = rna_input.rename_axis('sra')\
        .rename_axis(None, axis='columns').reset_index()

    # Factor names and Orthogroup names
    ortho_names = (rna_input.columns[1:]).tolist()

    # Merge (inner join) factors and RNAseq data by "sra".
    df = pd.merge(factor_input, rna_input, on='sra', how='inner')

    # Print out the after-merge shapes of dataframes
    print('Dataframe shape after merge: {}'.format(df.shape))

    # Transform and standardize RNAseq data (z-score by row)
    if logscale:
        df.iloc[:, fnc:] = (df.iloc[:, fnc:]).apply(lambda x: np.log2(x+1.0))

    if zscore:
        df.iloc[:, fnc:] = scale(df.iloc[:, fnc:], axis=1)

    if group:
        df_grouped = df.groupby('sample_id')
        df = pd.concat([df_grouped[factor_list].first(),
                        df_grouped[ortho_names].mean()], axis=1)

        # Print out the shapes of dataframe after grouping
        print('Data shape after grouping: {}'.format(df.shape))

    return df, ortho_names


def colorscale_from_matplotlib_cmap(cmap : object) -> list:
    """Generate colorscale for keppler mapper from a matplotlib cmap.

    Parameters
    ----------
    cmap : matplotlib cmap
        For example matplotlib.pyplot.cm.plasma

    Returns
    -------
    colorscale : list of tuples
        range of values from 0 to 1 inclusive zipped with RGB strings.
        colorscale format compatible with keppler mapper.

    """
    cmap_list = [cmap(el) for el in np.arange(cmap.N)]
    rgb_strings = ["rgb({}, {}, {})".format(
            int(255 * el[0]), int(255 * el[1]), int(255 * el[2])
        ) for el in cmap_list
    ]

    idx_scale = np.linspace(0, 1., num=cmap.N)

    return list(zip(idx_scale, rgb_strings))


def get_color_dict(df : object, factors : list, colorfile : str) -> dict:
    """Generate color mapping for factor levels

    Parameters
    ----------
    df : pandas dataframe
        dataframe that contains the factors (expected categorical) as columns.
    factors : list
        names of factors (expected categorical) to generate color mappings for.
    colorfile : str
        path to saved pickle file containing color_dict.

    Returns
    -------
    color_dict : dictionary
        dictionary mapping factor leveles to colors (both hex and RGB codes)

    """
    if os.path.isfile(colorfile):
        with open(colorfile, "rb") as fp:
            color_dict = pickle.load(fp)
    
    else:
        cmap = plt.get_cmap("tab20")
        color_dict = dict()
        for f in factors:
            color_dict[f] = dict()
            cols, labels = pd.factorize(df[f], sort=True)
            labels = labels.tolist()
            colors = []
            for i, l in enumerate(labels):
                c_hex = to_hex(cmap(i))
                c_rgb = list(to_rgb(cmap(i)))
                color_dict[f][l] = {"hex": c_hex, "rgb": c_rgb}

        with open(colorfile, "wb") as fp:
            pickle.dump(obj=color_dict, file=fp, protocol=pickle.HIGHEST_PROTOCOL)
            print("wrote file to: {}".format(colorfile))

    return color_dict


def create_pie_chart_json(graph : object, df : object, factors : list,
                          filterfactor : str, filterlevel : str,
                          color_dict : dict, json_dir : str):
    """Generate a json file for mapper with pie-chart nodes

    Parameters
    ----------
    graph: dict
        graph dictionary returned by map method of kepler mapper
    df : pandas dataframe
        dataframe that contains the factors (expected categorical) as columns.
    factors : list
        names of factors (expected categorical) to generate color mappings for.
    filterfactor : str
        factor used to construct the lens ("tissue" or "stress" or "family")
    filterlevel : str
        factor level used to create idealized linear model.
        eg. "healthy" if filterfactor is "stress", "leaf" or "root" if
        filterfactor is "tissue".
    color_dict : dict
        dictionary containing factor levels to colors mapping
    json_dir : str
        path to the dictionary where json files are to be saved

    Returns
    -------
    None

    """
    lens = np.asarray(graph["lens"])
    for f in factors:
        vals, unique = pd.factorize(df[f], sort=True)
        json_dict = {"nodes": [], "links": []}
        json_dict["labels"] = list(color_dict[f].keys())
        json_dict["colors"] = [color_dict[f][l]["hex"] 
                                for l in json_dict["labels"]]
        json_dict["label_idx"] = list(range(len(json_dict["labels"])))

        for i, node_ID in enumerate(graph["nodes"]):
            testlist = vals[graph['nodes'][node_ID]].tolist()
            n = {
                "id": node_ID, # this is used for links
                "group": int(stats.mode(testlist)[0][0]),
                "node_size": len(vals[graph['nodes'][node_ID]]),
                "lens_mean": np.mean(lens[graph['nodes'][node_ID]].tolist()),
                "pieChart": []
                }
                
            for j in set(testlist):
                slice = {"index": j,
                        "color": json_dict["colors"][j],
                        "percent": 100*testlist.count(j)/len(testlist)}
                n["pieChart"].append(slice)

            json_dict["nodes"].append(n)

        # we need a double for loop since some link values are length > 1
        for k in range(0, len(graph['links'])):
            for l in list(graph['links'].values())[k]:
                source = list(graph['links'].keys())[k]
                target = l
                value = len(set(df.loc[graph['nodes'][source]]['sra'].tolist()) & set(df.loc[graph['nodes'][target]]['sra'].tolist()))
                json_dict["links"].append({"source": source, "target": target, "value": value})
        
        json_file = f"{f}_data_{filterfactor}_{filterlevel}_lens.json"
        json_path = json_dir + "/" + json_file
        with open(json_path, 'w') as outfile:
            json.dump(json_dict, outfile)

        print(json_dict["labels"])
        print(json_dict["colors"])
    
    return None


def create_grouped_pie_chart_json(graph : object, df : object, 
                                  gdata : object, factors : list,
                                  filterfactor : str, filterlevel : str,
                                  color_dict : dict, json_dir : str):
    """Generate a json file for grouped mapper with pie-chart nodes

    Parameters
    ----------
    graph: dict
        graph dictionary returned by map method of kepler mapper
    df : pandas dataframe
        dataframe that contains the factors (expected categorical) as columns.
    gdata : pandas dataframe
        dataframe that contains group assignments for mapper nodes.
    factors : list
        names of factors (expected categorical) to generate color mappings for.
    filterfactor : str
        factor used to construct the lens ("tissue" or "stress" or "family")
    filterlevel : str
        factor level used to create idealized linear model.
        eg. "healthy" if filterfactor is "stress", "leaf" or "root" if
        filterfactor is "tissue".
    color_dict : dict
        dictionary containing factor levels to colors mapping
    json_dir : str
        path to the dictionary where json files are to be saved

    Returns
    -------
    None

    """

    lens = graph["lens"]
    gdf = gdata.groupby(["GroupIndex", "Group"]).size().reset_index()
    for f in factors:
        vals, unique = pd.factorize(df[f], sort=True)
        lens = np.asarray(graph["lens"]).reshape(-1, 1)
        json_dict = {"nodes": [], "links": []}
        json_dict["labels"] = unique.tolist()
        json_dict["label_idx"] = list(range(len(json_dict["labels"])))
        json_dict["colors"] = [color_dict[f][l]["hex"] for l in json_dict["labels"]]

        for tup in gdf.itertuples():
            idx = tup.GroupIndex
            name = tup.Group
            clusters = gdata[gdata["GroupIndex"]==idx]["NodeID"].tolist()
            samples = []
            for c in clusters:
                samples += graph["nodes"][c]

            samples = list(set(samples))
            testlist = vals[samples].tolist()
            n = {
                "id": idx,
                "name": name,
                "group": int(stats.mode(testlist)[0][0]),
                "node_size": len(samples),
                "samples": samples,
                "lens_mean": np.mean(lens[samples].tolist()),
                "pieChart": []
                }
                
            for j in set(testlist):
                slice = {"index": j,
                        "color": json_dict["colors"][j],
                        "percent": 100*testlist.count(j)/len(testlist)}
                n["pieChart"].append(slice)

            json_dict["nodes"].append(n)

        nodeIDs = gdf["GroupIndex"].ravel()
        for pair in combinations(nodeIDs, 2):
            set1 = set(json_dict["nodes"][pair[0] - 1]["samples"])
            set2 = set(json_dict["nodes"][pair[1] - 1]["samples"])
            linkset = set1.intersection(set2)
            if len(linkset) > 0:
                link = {"source": int(pair[0]),
                        "target": int(pair[1]),
                        "value": len(linkset),
                        "samples": list(linkset)
                    }
                json_dict["links"].append(link)

        json_file = f"grouped_{f}_data_{filterfactor}_{filterlevel}_lens.json"
        json_path = json_dir + "/" + json_file
        with open(json_path, 'w') as outfile:
            json.dump(json_dict, outfile)

        print("Wrote json file to {}".format(json_path))
    
    return None



if __name__ == "__main__":

    datadir = '../data/SVA_Corrected_20210928'
    factorfile = datadir + '/clean_metadata.csv'
    rnafile = datadir + '/clean_RNAseq_sv_corrected.csv'
    factors = ['stress', 'tissue', 'family']
    levels = ['healthy', 'leaf', 'Poaceae']

    data, orthos = loaddata(factorfile, rnafile, factors)
    print(data.shape)
