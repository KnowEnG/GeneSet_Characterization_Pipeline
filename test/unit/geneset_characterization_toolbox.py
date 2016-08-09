"""
Created on Tue Jun 28 14:39:35 2016
@author: The Gene Sets Characterization dev team
"""

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.preprocessing import normalize

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

def perform_fisher_exact_test(
        prop_gene_network_sparse, reverse_prop_gene_network_n1_names_dict,
        spreadsheet_df, results_dir):
    """ central loop: compute components for fisher exact test.
    Args:
        prop_gene_network_sparse: sparse matrix of network gene set.
        reverse_prop_gene_network_n1_names_dict: look up table of sparse matrix.
        spreadsheet_df: the dataframe of user gene set.
        results_dir: directory name to write results.
    """
    df_col = ["user gene", "property", "count", "user count", "gene count", "overlap", "pval"]
    gene_count = prop_gene_network_sparse.sum(axis=0)
    universe_count = spreadsheet_df.shape[0]
    df_val = []

    col_list = spreadsheet_df.columns.values
    for col in col_list:
        new_user_set = spreadsheet_df.loc[:, col]
        user_count = np.sum(new_user_set.values)
        overlap_count = prop_gene_network_sparse.T.dot(new_user_set.values)

        for i in range(0, len(reverse_prop_gene_network_n1_names_dict)):
            table = build_fisher_contigency_table(
                overlap_count[i], user_count, gene_count[0, i], universe_count)
            oddsratio, pvalue = stats.fisher_exact(table, alternative="greater")

            if overlap_count[i] != 0:
                row_item = [col, reverse_prop_gene_network_n1_names_dict[i], int(universe_count),
                            int(user_count), int(gene_count[0, i]), int(overlap_count[i]), round(pvalue, 3)]
                df_val.append(row_item)

    result_df = pd.DataFrame(df_val, columns=df_col).sort_values("pval", ascending=1)
    return result_df



def perform_DRaWR(network_sparse, spreadsheet_df, len_gene_names, run_parameters):
    """ calculate random walk with global network and user set gene sets  and write output.
    Args:
        network_sparse: sparse matrix of global network.
        spreadsheet_df: dataframe of user gene sets.
        len_gene_names: length of genes in the in the user spreadsheet.
        run_parameters: parameters dictionary.
    """
    hetero_network = normalize(network_sparse, norm='l1', axis=0)
    new_spreadsheet_df = append_column_to_spreadsheet(spreadsheet_df, len_gene_names)

    final_spreadsheet_matrix, step = smooth_matrix_with_rwr(
        normalize(new_spreadsheet_df, norm='l1', axis=0), hetero_network, run_parameters)

    final_spreadsheet_df = pd.DataFrame(
        final_spreadsheet_matrix, index=new_spreadsheet_df.index.values,
        columns=new_spreadsheet_df.columns.values)

    final_spreadsheet_df = final_spreadsheet_df.iloc[len_gene_names:]
    for col in final_spreadsheet_df.columns.values[:-1]:
        final_spreadsheet_df[col] = final_spreadsheet_df[col] - final_spreadsheet_df['base']
        final_spreadsheet_df[col] = final_spreadsheet_df.sort_values(col, ascending=0).index.values

    final_spreadsheet_df['base'] = \
        final_spreadsheet_df.sort_values('base', ascending=0).index.values

    return

def run_fisher(run_parameters):
    ''' wrapper: call sequence to perform fisher gene-set characterization
    Args:
        run_parameters: dictionary of run parameters
    '''
    # -----------------------------------
    # - Data read and extraction Section -
    # -----------------------------------
    spreadsheet_df = get_spreadsheet_df(run_parameters)
    prop_gene_network_df = get_network_df(run_parameters['pg_network_file_name'])

    spreadsheet_gene_names = extract_spreadsheet_gene_names(spreadsheet_df)

    prop_gene_network_n1_names,\
    prop_gene_network_n2_names = extract_network_node_names(prop_gene_network_df)

    # -----------------------------------------------------------------------
    # - limit the gene set to the intersection of network and user gene set -
    # -----------------------------------------------------------------------
    common_gene_names = find_common_node_names(prop_gene_network_n2_names, spreadsheet_gene_names)

    common_gene_names_dict = create_node_names_dict(common_gene_names)

    prop_gene_network_n1_names_dict = create_node_names_dict(prop_gene_network_n1_names)

    reverse_prop_gene_network_n1_names_dict = create_reverse_node_names_dict(
        prop_gene_network_n1_names_dict)

    # ----------------------------------------------------------------------------
    # - restrict spreadsheet and network to common genes and drop everthing else -
    # ----------------------------------------------------------------------------
    spreadsheet_df = update_spreadsheet_df(spreadsheet_df, common_gene_names)
    prop_gene_network_df = update_network_df(prop_gene_network_df, common_gene_names, "node_2")

    # ----------------------------------------------------------------------------
    # - map every gene name to an integer index in sequential order startng at 0 -
    # ----------------------------------------------------------------------------
    prop_gene_network_df = map_node_names_to_index(
        prop_gene_network_df, prop_gene_network_n1_names_dict, "node_1")
    prop_gene_network_df = map_node_names_to_index(
        prop_gene_network_df, common_gene_names_dict, "node_2")

    # --------------------------------------------
    # - store the network in a csr sparse format -
    # --------------------------------------------
    universe_count = len(common_gene_names)
    prop_gene_network_sparse = convert_network_df_to_sparse(
        prop_gene_network_df, universe_count, len(prop_gene_network_n1_names))
    perform_fisher_exact_test(
        prop_gene_network_sparse, reverse_prop_gene_network_n1_names_dict,
        spreadsheet_df, run_parameters['results_directory'])

    return
def run_DRaWR(run_parameters):
    ''' wrapper: call sequence to perform random walk with restart
    Args:
        run_parameters: dictionary of run parameters
    '''
    spreadsheet_df = get_spreadsheet_df(run_parameters)
    pg_network_df = get_network_df(run_parameters['pg_network_file_name'])
    gg_network_df = get_network_df(run_parameters['gg_network_file_name'])

    pg_network_n1_names,\
    pg_network_n2_names = extract_network_node_names(pg_network_df)

    gg_network_n1_names,\
    gg_network_n2_names = extract_network_node_names(gg_network_df)

    # limit the gene set to the intersection of networks (gene_gene and prop_gene) and user gene set
    unique_gene_names = find_unique_node_names(gg_network_n1_names, gg_network_n2_names)
    unique_gene_names = find_unique_node_names(unique_gene_names, pg_network_n2_names)
    unique_all_node_names = unique_gene_names + pg_network_n1_names
    unique_gene_names_dict = create_node_names_dict(unique_gene_names)
    pg_network_n1_names_dict = create_node_names_dict(
        pg_network_n1_names, len(unique_gene_names))

    # restrict spreadsheet to unique genes and drop everthing else
    spreadsheet_df = update_spreadsheet_df(spreadsheet_df, unique_all_node_names)
    # map every gene name to a sequential integer index
    gg_network_df = map_node_names_to_index(gg_network_df, unique_gene_names_dict, "node_1")
    gg_network_df = map_node_names_to_index(gg_network_df, unique_gene_names_dict, "node_2")
    pg_network_df = map_node_names_to_index(pg_network_df, pg_network_n1_names_dict, "node_1")
    pg_network_df = map_node_names_to_index(pg_network_df, unique_gene_names_dict, "node_2")

    gg_network_df = symmetrize_df(gg_network_df)
    pg_network_df = symmetrize_df(pg_network_df)

    gg_network_df = normalize_df(gg_network_df, 'wt')
    pg_network_df = normalize_df(pg_network_df, 'wt')

    hybrid_network_df = form_hybrid_network_df([gg_network_df, pg_network_df])

    # store the network in a csr sparse format
    network_sparse = convert_network_df_to_sparse(
        hybrid_network_df, len(unique_all_node_names), len(unique_all_node_names))

    perform_DRaWR(network_sparse, spreadsheet_df, len(unique_gene_names), run_parameters)

    return
