"""
Created on Tue Jun 28 14:39:35 2016
@author: The Gene Sets Characterization dev team
"""
import os
import time
import numpy as np
import numpy.linalg as LA
import pandas as pd
from scipy import linalg
from scipy import stats
from sklearn.preprocessing import normalize
from sklearn.metrics.pairwise import cosine_similarity
import knpackage.toolbox as kn

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
    reverse_prop_dict = kn.create_reverse_node_names_dict(prop_gene_network_n1_names_dict)
    # ----------------------------------------------------------------------------
    # - restrict spreadsheet and network to common genes and drop everthing else -
    # ----------------------------------------------------------------------------
    new_spreadsheet_df = kn.update_spreadsheet_df(spreadsheet_df, common_gene_names)
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
    fisher_contingency_pval = get_fisher_exact_test(prop_gene_network_sparse,
        reverse_prop_dict, new_spreadsheet_df)

    save_fisher_test_result(fisher_contingency_pval, 
        run_parameters['results_directory'], spreadsheet_df.columns.values)
    map_and_save_droplist(spreadsheet_df, common_gene_names, 'fisher_droplist', run_parameters)

    return fisher_contingency_pval

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

    new_spreadsheet_df = kn.update_spreadsheet_df(spreadsheet_df, unique_all_node_names)
    
    base_col = np.append(np.ones(len(unique_gene_names), dtype=np.int),
                         np.zeros(len(set(pg_network_n1_names)), dtype=np.int))

    new_spreadsheet_df = kn.append_column_to_spreadsheet(new_spreadsheet_df, base_col, 'base')

    final_spreadsheet_df = get_DRaWR(network_sparse, new_spreadsheet_df, 
        len(unique_gene_names), run_parameters)

    save_timestamped_df(final_spreadsheet_df, run_parameters['results_directory'], 'DRaWR_result')
    map_and_save_droplist(spreadsheet_df, unique_gene_names, 'DRaWR_droplist', run_parameters)
    return final_spreadsheet_df

def run_net_path(run_parameters):
    ''' wrapper: call sequence to perform net path
    Args:
        run_parameters: dictionary of run parameters
    '''
    spreadsheet_df = kn.get_spreadsheet_df(run_parameters['spreadsheet_name_full_path'])
    pg_network_df = kn.get_network_df(run_parameters['pg_network_name_full_path'])
    gg_network_df = kn.get_network_df(run_parameters['gg_network_name_full_path'])

    pg_network_n1_names, \
    pg_network_n2_names = kn.extract_network_node_names(pg_network_df)
    gg_network_n1_names, \
    gg_network_n2_names = kn.extract_network_node_names(gg_network_df)

    # limit the gene set to the intersection of networks (gene_gene and prop_gene) and user gene set
    unique_gene_names = kn.find_unique_node_names(gg_network_n1_names, gg_network_n2_names)
    unique_all_node_names = unique_gene_names + pg_network_n1_names
    pg_network_df = kn.update_network_df(pg_network_df, unique_gene_names, 'node_2')

    unique_gene_names_dict = kn.create_node_names_dict(unique_gene_names)
    pg_network_n1_names_dict = kn.create_node_names_dict(
        pg_network_n1_names, len(unique_gene_names))


    spreadsheet_df = kn.update_spreadsheet_df(spreadsheet_df, unique_gene_names)
    gg_network_df = kn.map_node_names_to_index(gg_network_df, unique_gene_names_dict, "node_1")
    gg_network_df = kn.map_node_names_to_index(gg_network_df, unique_gene_names_dict, "node_2")
    pg_network_df = kn.map_node_names_to_index(pg_network_df, pg_network_n1_names_dict, "node_1")
    pg_network_df = kn.map_node_names_to_index(pg_network_df, unique_gene_names_dict, "node_2")

    gg_network_df = kn.symmetrize_df(gg_network_df)
    pg_network_df = kn.symmetrize_df(pg_network_df)

    hybrid_network_df = kn.form_hybrid_network_df([gg_network_df, pg_network_df])
    network_sparse = kn.convert_network_df_to_sparse(
        hybrid_network_df, len(unique_all_node_names), len(unique_all_node_names))

    ret = perform_net_path(
        spreadsheet_df, network_sparse, unique_gene_names, pg_network_n1_names, run_parameters)

    save_timestamped_df(ret, run_parameters['results_directory'], 'net_path_result')
    map_and_save_droplist(spreadsheet_df, unique_gene_names, 'net_path_droplist', run_parameters)
    
    return ret


def perform_k_SVD(smooth_spreadsheet_matrix, k):
    """Perform SVD on input matrix.

    Args:
        smooth_spreadsheet_matrix: Input matrix to perform SVD on.
        k: Number of singular values and vectors to compute.

    Returns:
        U_unitary_matrix: Unitary matrix having left singular vectors as column.
        S_full_squared_matrix: Matrix with diagonal to be k singular values.
    """
    U_unitary_matrix, singular_value, V_unitary_matrix = linalg.svd(smooth_spreadsheet_matrix)
    S_full_squared_matrix = np.zeros((k, k))
    np.fill_diagonal(S_full_squared_matrix, np.sqrt(singular_value[:k]))
    U_unitary_matrix = U_unitary_matrix[:, :k]
    return U_unitary_matrix, S_full_squared_matrix

def project_matrix_to_new_space_and_split(U_unitary_matrix, S_full_squared_matrix,
                                          unique_gene_length):
    """This to project matrix to the new space and split into gene and
    property two parts.

    Args:
        U_unitary_matrix: Unitary matrix having left singular vectors as column.
        S_full_squared_matrix: Matrix with diagonal to be k singular values.
        unique_gene_length: Length of unique genes.

    Returns:
        g_newspace_matrix: gene matrix projected to the new space.
        p_newspace_matrix: property matrix projected to new space.
    """
    L = U_unitary_matrix.dot(S_full_squared_matrix)
    g_newspace_matrix = L[:unique_gene_length]
    p_newspace_matrix = L[unique_gene_length:]

    return g_newspace_matrix, p_newspace_matrix

def perform_cosine_correlation(g_newspace_matrix, p_newspace_matrix,
                               gene_names, property_name):
    """This is to perform cosine similarity on two input matrixes.

    Args:
        g_newspace_matrix: input gene matrix.
        p_newspace_matrix: intput property matrix.
        gene_names: list of gene names.
        property_name: list of property names.

    Returns:
        cosine_matrix_df: dataframe with cosine values calculated from input
        matrixes.
    """
    cosine_matrix = cosine_similarity(g_newspace_matrix, p_newspace_matrix)
    cosine_matrix_df = pd.DataFrame(cosine_matrix, index=gene_names, columns=property_name)
    return cosine_matrix_df

def smooth_final_spreadsheet_matrix(final_rwr_matrix, gene_length):
    """This is to add pseudo count to the input matrix.

    Args:
        final_rwr_matrix: input matrix.
        gene_length: length of genes.

    Returns:
        smooth_rwr_matrix: the smoothed matrix with pseudo counts.
    """
    assert ((final_rwr_matrix >= 0).all())
    eps = np.float(1/gene_length)
    smooth_rwr_matrix = np.log(final_rwr_matrix + eps) - np.log(eps)
    return smooth_rwr_matrix

def rank_property(spreadsheet_df, cosine_matrix_df):
    """This is to rank property based on cosine values:

    Args:
        spreadsheet_df: user supplied spreadsheet dataframe.
        cosine_matrix_df: dataframe with cosine value.

    Returns:
        property_rank_df: dataframe with ranked property names in each column.
    """
    property_rank_df = pd.DataFrame(columns=spreadsheet_df.columns.values)
    for col_name in spreadsheet_df.columns.values:
        user_gene_list = spreadsheet_df[spreadsheet_df[col_name] == 1].index.values
        new_spreadsheet_df = cosine_matrix_df.loc[user_gene_list].sum()
        property_rank_df[col_name] = new_spreadsheet_df.sort_values(ascending=False).index.values
    return property_rank_df

def save_cosine_matrix_df(cosine_matrix_df, run_parameters):
    """This is to save the cosine matrix df to output file

    Args:
        cosine_matrix_df: dataframe with cosine value.
        run_parameters: parameters dictionary.
    """
    new_file_name = kn.create_timestamped_filename("cosine_matrix", "df")
    cosine_matrix_df.to_csv(os.path.join(run_parameters['results_directory'], new_file_name), header=True, index=True, sep='\t')

def perform_net_path(spreadsheet_df, network_sparse, unique_gene_names,
                     pg_network_n1_names, run_parameters):
    """Perform net path method on gene characterization.

    Args:
        network_sparse: sparse matrix of global network.
        spreadsheet_df: dataframe of user gene sets.
        unique_gene_names: list of genes in the in the user spreadsheet.
        pg_network_n1_names: list of property names in the network.
        run_parameters: parameters dictionary.

    Returns:
        property_rank_df: dataframe with ranked property names in each column.
    """
    hetero_network = normalize(network_sparse, norm='l1', axis=0)
    restart = np.eye(hetero_network.shape[0])
    final_rwr_matrix, step = kn.smooth_matrix_with_rwr(
        restart, hetero_network, run_parameters)
    smooth_rwr_matrix = smooth_final_spreadsheet_matrix(final_rwr_matrix, len(unique_gene_names))
    unique_all_node_names = unique_gene_names + pg_network_n1_names
    U_unitary_matrix, S_full_squared_matrix = perform_k_SVD(smooth_rwr_matrix, 
        int(run_parameters['k_space']), unique_all_node_names)
    g_newspace_matrix, p_newspace_matrix = project_matrix_to_new_space_and_split(
        U_unitary_matrix, S_full_squared_matrix, len(unique_gene_names))
    cosine_matrix_df = perform_cosine_correlation(
        g_newspace_matrix, p_newspace_matrix, unique_gene_names, pg_network_n1_names)
    save_cosine_matrix_df(cosine_matrix_df, run_parameters)

    property_rank_df = rank_property(spreadsheet_df, cosine_matrix_df)
    return property_rank_df

def build_fisher_contingency_table(overlap_count, user_count, gene_count, count):
    """ build contingency table for fisher exact test.

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

def get_fisher_exact_test(prop_gene_network_sparse, sparse_dict, spreadsheet_df):
    """ central loop: compute components for fisher exact test.

    Args:
        prop_gene_network_sparse: sparse matrix of network gene set.
        sparse_dict: look up table of sparse matrix.
        spreadsheet_df: the dataframe of user gene set.

    Returns:
        fisher_contingency_pval: list of seven items lists.
    """
    universe_count = spreadsheet_df.shape[0]
    overlap_count = prop_gene_network_sparse.T.dot(spreadsheet_df.values)
    user_count = np.sum(spreadsheet_df.values, axis=0)
    gene_count = prop_gene_network_sparse.sum(axis=0)
    set_list = spreadsheet_df.columns.values
    fisher_contingency_pval = []

    for i in range(overlap_count.shape[0]):
        for j in range(overlap_count.shape[1]):
            table = build_fisher_contingency_table(overlap_count[i, j], user_count[j],
                                                  gene_count[0, i], universe_count)
            pvalue = stats.fisher_exact(table, alternative="greater")[1]
            row_item = [set_list[j], sparse_dict[i], int(universe_count), int(user_count[j]),
                        int(gene_count[0, i]), int(overlap_count[i, j]), pvalue]
            fisher_contingency_pval.append(row_item)

    return fisher_contingency_pval

def get_DRaWR(network_sparse, new_spreadsheet_df, len_gene_names, run_parameters):
    """ calculate random walk with global network and user set gene sets  and write output.
    Args:
        network_sparse: sparse matrix of global network.
        new_spreadsheet_df: dataframe of user gene sets.
        len_gene_names: length of genes in the in the user spreadsheet.
        run_parameters: parameters dictionary.
    Returns:
        final_spreadsheet_df: dataframe with ranked property names.
    """
    hetero_network = normalize(network_sparse, norm='l1', axis=0)
    final_spreadsheet_matrix, step = kn.smooth_matrix_with_rwr(
        normalize(new_spreadsheet_df, norm='l1', axis=0), hetero_network, run_parameters)

    final_spreadsheet_df = pd.DataFrame(final_spreadsheet_matrix[len_gene_names:])
    final_spreadsheet_df.index = new_spreadsheet_df.index.values[len_gene_names:]
    final_spreadsheet_df.columns = new_spreadsheet_df.columns.values

    final_spreadsheet_df.iloc[:, :-1] = final_spreadsheet_df.iloc[:, :-1].apply(
        lambda x: (x - final_spreadsheet_df['base']).sort_values(ascending=0).index.values)
    final_spreadsheet_df = final_spreadsheet_df.drop('base', 1)

    return final_spreadsheet_df

def save_fisher_test_result(fisher_contingency_pval, results_dir, set_list):
    """ Save two output files of fisher exact test results.

    Args:
        fisher_contingency_pval: list of seven items lists.
        set_list: column values of spreadsheet.
    """
    df_col = ["user gene", "property", "count", "user count", "gene count", "overlap", "pval"]
    result_df = pd.DataFrame(fisher_contingency_pval, columns=df_col).sort_values("pval", ascending=1) 
    save_timestamped_df(result_df, results_dir, 'fisher_result')

    new_result_df = pd.DataFrame(columns=set_list)
    for gene_set in set_list:
        new_result_df.loc[:, gene_set] = result_df[result_df['user gene']==gene_set].values[:, 1] 

    save_timestamped_df(new_result_df, results_dir, 'fisher_result_geneset_property')

def save_timestamped_df(input_df, results_dir, name):
    """ Save dataframe to files with timestamped name.

    Args:
        fisher_contingency_pval: list of seven items lists.
        results_dir: directory to save outputs.
        name: file name.
    """
    file_name = kn.create_timestamped_filename(name, "df")
    kn.save_df(input_df, results_dir, file_name)

def map_and_save_droplist(spreadsheet_df, gene_names, droplist_name, run_parameters):
    """This is to map and save droplist

    Args:
        spreadsheet_df: user supplied spreadsheet dataframe.
        gene_names: list of genes.
        droplist_name: name of droplist file.
        run_parameters: dictionary of run parameters.

    Returns:
        property_rank_df: dataframe with ranked property names in each column.
    """
    droplist = kn.find_dropped_node_names(spreadsheet_df, gene_names)
    file_name = kn.create_timestamped_filename(droplist_name, "tsv")
    kn.save_df(pd.DataFrame(droplist, columns=[droplist_name]),
               run_parameters['results_directory'], file_name)
