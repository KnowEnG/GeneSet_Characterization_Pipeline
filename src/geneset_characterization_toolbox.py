"""
Created on Tue Jun 28 14:39:35 2016
@author: The Gene Sets Characterization dev team
"""

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.preprocessing import normalize
import knpackage.toolbox as kn

def build_fisher_contigency_table(overlap_count, user_count, gene_count, count):
    """ build contigency table for fisher exact test.

    Args:
        overlap_count: count of overlaps in user gene set and network gene set.
        user_count: count of ones in user gene set.
        gene_count: count of ones in network gene set
        count: number of universe genes.
    Returns:
        table: the contigency table used in fisher test.
    """
    table = np.zeros(shape=(2, 2))
    table[0, 0] = overlap_count
    table[0, 1] = user_count - overlap_count
    table[1, 0] = gene_count - overlap_count
    table[1, 1] = count - user_count - gene_count + overlap_count

    return table

def perform_fisher_exact_test(prop_gene_network_sparse, sparse_dict,
                              spreadsheet_df, results_dir):
    """ central loop: compute components for fisher exact test.

    Args:
        prop_gene_network_sparse: sparse matrix of network gene set.
        sparse_dict: look up table of sparse matrix.
        spreadsheet_df: the dataframe of user gene set.
        results_dir: directory name to write results.
    """
    universe_count = spreadsheet_df.shape[0]
    overlap_count = prop_gene_network_sparse.T.dot(spreadsheet_df.values)
    user_count = np.sum(spreadsheet_df.values, axis=0)
    gene_count = prop_gene_network_sparse.sum(axis=0)
    set_list = spreadsheet_df.columns.values
    df_val = []
    
    for i in range(overlap_count.shape[0]):
        for j in range(overlap_count.shape[1]):
            table = build_fisher_contigency_table(overlap_count[i, j], user_count[j],
                                                  gene_count[0, i], universe_count)
            pvalue = stats.fisher_exact(table, alternative="greater")[1]
            if overlap_count[i, j] != 0:
                row_item = [set_list[j], sparse_dict[i], int(universe_count), int(user_count[j]),
                            int(gene_count[0, i]), int(overlap_count[i, j]), pvalue]
                df_val.append(row_item)
    df_col = ["user gene", "property", "count", "user count", "gene count", "overlap", "pval"]
    result_df = pd.DataFrame(df_val, columns=df_col).sort_values("pval", ascending=1)
    file_name = kn.create_timestamped_filename("fisher_result", "df")
    kn.save_df(result_df, results_dir, file_name)

    return result_df

def perform_DRaWR(network_sparse, new_spreadsheet_df, len_gene_names, run_parameters):
    """ calculate random walk with global network and user set gene sets  and write output.
    Args:
        network_sparse: sparse matrix of global network.
        new_spreadsheet_df: dataframe of user gene sets.
        len_gene_names: length of genes in the in the user spreadsheet.
        run_parameters: parameters dictionary.
    """
    hetero_network = normalize(network_sparse, norm='l1', axis=0)
    final_spreadsheet_matrix, step = kn.smooth_matrix_with_rwr(
        normalize(new_spreadsheet_df, norm='l1', axis=0), hetero_network, run_parameters)

    final_spreadsheet_df = pd.DataFrame(final_spreadsheet_matrix[len_gene_names:])
    final_spreadsheet_df.index = new_spreadsheet_df.index.values[len_gene_names:]
    final_spreadsheet_df.columns = new_spreadsheet_df.columns.values

    final_spreadsheet_df.iloc[:, :-1] = final_spreadsheet_df.iloc[:, :-1].apply(
        lambda x: (x - final_spreadsheet_df['base']).sort_values(ascending=0).index.values)

    final_spreadsheet_df['base'] = \
        final_spreadsheet_df['base'].sort_values(ascending=0).index.values

    final_spreadsheet_df.index = range(final_spreadsheet_df.shape[0])
    file_name = kn.create_timestamped_filename("DRaWR_result", "df")
    kn.save_df(final_spreadsheet_df, run_parameters['results_directory'], file_name)

    return final_spreadsheet_df

def run_fisher(run_parameters):
    ''' wrapper: call sequence to perform fisher gene-set characterization
    Args:
        run_parameters: dictionary of run parameters
    '''
    # -----------------------------------
    # - Data read and extraction Section -
    # -----------------------------------
    spreadsheet_df = kn.get_spreadsheet_df(run_parameters['spreadsheet_name_full_path'])
    prop_gene_network_df = kn.get_network_df(run_parameters['pg_network_name_full_path'])

    spreadsheet_gene_names = kn.extract_spreadsheet_gene_names(spreadsheet_df)

    prop_gene_network_n1_names,\
    prop_gene_network_n2_names = kn.extract_network_node_names(prop_gene_network_df)

    # -----------------------------------------------------------------------
    # - limit the gene set to the intersection of network and user gene set -
    # -----------------------------------------------------------------------
    common_gene_names = kn.find_common_node_names(prop_gene_network_n2_names, spreadsheet_gene_names)

    common_gene_names_dict = kn.create_node_names_dict(common_gene_names)

    prop_gene_network_n1_names_dict = kn.create_node_names_dict(prop_gene_network_n1_names)

    reverse_prop_gene_network_n1_names_dict = kn.create_reverse_node_names_dict(
        prop_gene_network_n1_names_dict)

    # ----------------------------------------------------------------------------
    # - restrict spreadsheet and network to common genes and drop everthing else -
    # ----------------------------------------------------------------------------
    droplist = kn.find_dropped_node_names(spreadsheet_df, common_gene_names)
    file_name = kn.create_timestamped_filename("fisher_droplist", "tsv")
    kn.save_df(pd.DataFrame(droplist, columns=['droplist']),
               run_parameters['results_directory'], file_name)
    spreadsheet_df = kn.update_spreadsheet_df(spreadsheet_df, common_gene_names)
    prop_gene_network_df = kn.update_network_df(prop_gene_network_df, common_gene_names, "node_2")

    # ----------------------------------------------------------------------------
    # - map every gene name to an integer index in sequential order startng at 0 -
    # ----------------------------------------------------------------------------
    prop_gene_network_df = kn.map_node_names_to_index(
        prop_gene_network_df, prop_gene_network_n1_names_dict, "node_1")
    prop_gene_network_df = kn.map_node_names_to_index(
        prop_gene_network_df, common_gene_names_dict, "node_2")

    # --------------------------------------------
    # - store the network in a csr sparse format -
    # --------------------------------------------
    universe_count = len(common_gene_names)
    prop_gene_network_sparse = kn.convert_network_df_to_sparse(
        prop_gene_network_df, universe_count, len(prop_gene_network_n1_names))
    ret = perform_fisher_exact_test(
        prop_gene_network_sparse, reverse_prop_gene_network_n1_names_dict,
        spreadsheet_df, run_parameters['results_directory'])

    return ret
def run_DRaWR(run_parameters):
    ''' wrapper: call sequence to perform random walk with restart
    Args:
        run_parameters: dictionary of run parameters
    '''
    spreadsheet_df = kn.get_spreadsheet_df(run_parameters['spreadsheet_name_full_path'])
    pg_network_df = kn.get_network_df(run_parameters['pg_network_name_full_path'])
    gg_network_df = kn.get_network_df(run_parameters['gg_network_name_full_path'])

    pg_network_n1_names,\
    pg_network_n2_names = kn.extract_network_node_names(pg_network_df)

    gg_network_n1_names,\
    gg_network_n2_names = kn.extract_network_node_names(gg_network_df)

    # limit the gene set to the intersection of networks (gene_gene and prop_gene) and user gene set
    unique_gene_names = kn.find_unique_node_names(gg_network_n1_names, gg_network_n2_names)
    unique_gene_names = kn.find_unique_node_names(unique_gene_names, pg_network_n2_names)
    unique_all_node_names = unique_gene_names + pg_network_n1_names
    unique_gene_names_dict = kn.create_node_names_dict(unique_gene_names)
    pg_network_n1_names_dict = kn.create_node_names_dict(
        pg_network_n1_names, len(unique_gene_names))

    # restrict spreadsheet to unique genes and drop everthing else
    droplist = kn.find_dropped_node_names(spreadsheet_df, unique_gene_names)
    file_name = kn.create_timestamped_filename("DRaWR_droplist", "tsv")
    kn.save_df(pd.DataFrame(droplist, columns=['droplist']),
               run_parameters['results_directory'], file_name)
    spreadsheet_df = kn.update_spreadsheet_df(spreadsheet_df, unique_all_node_names)
    # map every gene name to a sequential integer index
    gg_network_df = kn.map_node_names_to_index(gg_network_df, unique_gene_names_dict, "node_1")
    gg_network_df = kn.map_node_names_to_index(gg_network_df, unique_gene_names_dict, "node_2")
    pg_network_df = kn.map_node_names_to_index(pg_network_df, pg_network_n1_names_dict, "node_1")
    pg_network_df = kn.map_node_names_to_index(pg_network_df, unique_gene_names_dict, "node_2")

    gg_network_df = kn.symmetrize_df(gg_network_df)
    pg_network_df = kn.symmetrize_df(pg_network_df)

    gg_network_df = kn.normalize_network_df_by_sum(gg_network_df, 'wt')
    pg_network_df = kn.normalize_network_df_by_sum(pg_network_df, 'wt')

    hybrid_network_df = kn.form_hybrid_network_df([gg_network_df, pg_network_df])

    # store the network in a csr sparse format
    network_sparse = kn.convert_network_df_to_sparse(
        hybrid_network_df, len(unique_all_node_names), len(unique_all_node_names))

    property_size = spreadsheet_df.shape[0] - len(unique_gene_names)
    base_col = np.append(np.ones(len(unique_gene_names), dtype=np.int),
                         np.zeros(property_size, dtype=np.int))

    new_spreadsheet_df = kn.append_column_to_spreadsheet(spreadsheet_df, base_col, 'base')

    ret = perform_DRaWR(network_sparse, new_spreadsheet_df, len(unique_gene_names), run_parameters)

    return ret
