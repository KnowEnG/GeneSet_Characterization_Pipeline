"""
Created on Tue Jun 28 14:39:35 2016
@author: The Gene Sets Characterization dev team
"""
import os
import sys
import numpy as np
import pandas as pd
from scipy import linalg
from scipy import stats
from sklearn.preprocessing import normalize
from sklearn.metrics.pairwise import cosine_similarity
import knpackage.toolbox as kn
import multiprocessing
import itertools
import knpackage.distributed_computing_utils as dstutil


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

    prop_gene_network_n1_names, \
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
    prop_gene_network_df['wt'] = 1
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
    fisher_contingency_pval = get_fisher_exact_test(
        prop_gene_network_sparse, reverse_prop_dict, new_spreadsheet_df)
    fisher_final_result = save_fisher_test_result(
        fisher_contingency_pval, run_parameters['results_directory'], spreadsheet_df.columns.values, 2)
    map_and_save_droplist(spreadsheet_df, common_gene_names, 'fisher_droplist', run_parameters)

    return fisher_final_result


def run_DRaWR(run_parameters):
    ''' wrapper: call sequence to perform random walk with restart
    Args:
        run_parameters: dictionary of run parameters
    '''

    network_sparse, unique_gene_names, \
    pg_network_n1_names = build_hybrid_sparse_matrix(run_parameters, True, True)

    unique_all_node_names = unique_gene_names + pg_network_n1_names
    spreadsheet_df = kn.get_spreadsheet_df(run_parameters['spreadsheet_name_full_path'])
    new_spreadsheet_df = kn.update_spreadsheet_df(spreadsheet_df, unique_all_node_names)

    unique_genes_length = len(unique_gene_names)
    property_length = len(set(pg_network_n1_names))
    base_col = np.append(np.ones(unique_genes_length, dtype=np.int),
                         np.zeros(property_length, dtype=np.int))

    new_spreadsheet_df = kn.append_column_to_spreadsheet(new_spreadsheet_df, base_col, 'base')
    hetero_network = normalize(network_sparse, norm='l1', axis=0)

    final_spreadsheet_matrix, step = kn.smooth_matrix_with_rwr(
        normalize(new_spreadsheet_df, norm='l1', axis=0), hetero_network, run_parameters)

    final_spreadsheet_df = pd.DataFrame(final_spreadsheet_matrix)
    final_spreadsheet_df.index = new_spreadsheet_df.index.values
    final_spreadsheet_df.columns = new_spreadsheet_df.columns.values
    prop_spreadsheet_df = rank_drawr_property(final_spreadsheet_df, pg_network_n1_names)

    spreadsheet_df_mask = final_spreadsheet_df.loc[final_spreadsheet_df.index.isin(spreadsheet_df.index)]
    gene_result_df = construct_drawr_result_df(
        spreadsheet_df_mask, 0, spreadsheet_df_mask.shape[0], True, run_parameters)
    prop_result_df = construct_drawr_result_df(
        final_spreadsheet_df, unique_genes_length, final_spreadsheet_df.shape[0], False, run_parameters)

    save_timestamped_df(prop_spreadsheet_df, run_parameters['results_directory'], 'DRaWR_ranked_by_property')
    save_timestamped_df(
        gene_result_df, run_parameters['results_directory'], 'DRaWR_sorted_by_gene_score')
    save_timestamped_df(
        prop_result_df, run_parameters['results_directory'], 'DRaWR_sorted_by_property_score')

    map_and_save_droplist(spreadsheet_df, unique_gene_names, 'DRaWR_droplist', run_parameters)

    return prop_spreadsheet_df


def run_net_path(run_parameters):
    ''' wrapper: call sequence to perform net path
    Args:
        run_parameters: dictionary of run parameters
    '''
    network_sparse, unique_gene_names, \
    pg_network_n1_names = build_hybrid_sparse_matrix(run_parameters, False, False)

    spreadsheet_df = kn.get_spreadsheet_df(run_parameters['spreadsheet_name_full_path'])
    new_spreadsheet_df = kn.update_spreadsheet_df(spreadsheet_df, unique_gene_names)

    hetero_network = normalize(network_sparse, norm='l1', axis=0)
    final_rwr_matrix, step = kn.smooth_matrix_with_rwr(
        np.eye(hetero_network.shape[0]), hetero_network, run_parameters)
    smooth_rwr_matrix = smooth_final_spreadsheet_matrix(final_rwr_matrix, len(unique_gene_names))

    cosine_matrix = get_net_path_results(len(unique_gene_names), smooth_rwr_matrix, run_parameters)

    cosine_matrix_df = pd.DataFrame(cosine_matrix, index=unique_gene_names, columns=pg_network_n1_names)
    # save_cosine_matrix_df(cosine_matrix_df, run_parameters)

    property_rank_df = rank_netpath_property(new_spreadsheet_df, cosine_matrix_df)
    prop_result_df = construct_netpath_result_df(new_spreadsheet_df, cosine_matrix_df)

    save_timestamped_df(property_rank_df, run_parameters['results_directory'], 'net_path_ranked_by_property')
    save_timestamped_df(
        prop_result_df, run_parameters['results_directory'], 'net_path_sorted_by_property_score')
    map_and_save_droplist(spreadsheet_df, unique_gene_names, 'net_path_droplist', run_parameters)

    return property_rank_df


def calculate_k_SVD(smooth_spreadsheet_matrix, k):
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


def smooth_final_spreadsheet_matrix(final_rwr_matrix, gene_length):
    """This is to add pseudo count to the input matrix.

    Args:
        final_rwr_matrix: input matrix.
        gene_length: length of genes.

    Returns:
        smooth_rwr_matrix: the smoothed matrix with pseudo counts.
    """
    assert ((final_rwr_matrix >= 0).all())
    eps = np.float(1 / gene_length)
    smooth_rwr_matrix = np.log(final_rwr_matrix + eps) - np.log(eps)
    return smooth_rwr_matrix


def save_cosine_matrix_df(cosine_matrix_df, run_parameters):
    """This is to save the cosine matrix df to output file

    Args:
        cosine_matrix_df: dataframe with cosine value.
        run_parameters: parameters dictionary.
    """
    new_file_name = kn.create_timestamped_filename("cosine_matrix", "df")
    cosine_matrix_df.to_csv(
        os.path.join(run_parameters['results_directory'], new_file_name), header=True, index=True, sep='\t')


def get_net_path_results(gene_length, smooth_rwr_matrix, run_parameters):
    """Perform net path method on gene characterization.

    Args:
        gene_length: length of gene list.
        smooth_rwr_matrix: smoothed matrix need to perform SVD on
        run_parameters: parameters dictionary.

    Returns:
        cosine_matrix: matrix with property and gene cosine similarity values.
    """

    U_unitary_matrix, S_full_squared_matrix = calculate_k_SVD(
        smooth_rwr_matrix, int(run_parameters['k_space']))

    g_newspace_matrix, p_newspace_matrix = project_matrix_to_new_space_and_split(
        U_unitary_matrix, S_full_squared_matrix, gene_length)

    cosine_matrix = cosine_similarity(g_newspace_matrix, p_newspace_matrix)

    return cosine_matrix


def rank_netpath_property(spreadsheet_df, cosine_matrix_df):
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


def construct_netpath_result_df(spreadsheet_df, cosine_matrix_df):
    """Construct a three-column netpath result dataframe with header as user_gene_set',
     'property_gene_set' and 'cosine_sum'

    Args:
        spreadsheet_df: user supplied spreadsheet dataframe.
        cosine_matrix_df: dataframe with cosine value.
    Returns:
        result_df: result five-column dataframe.
    """
    property_rank_df = pd.DataFrame(
        columns=spreadsheet_df.columns.values, index=cosine_matrix_df.columns.values)

    for col_name in spreadsheet_df.columns.values:
        user_gene_list = spreadsheet_df[spreadsheet_df[col_name] == 1].index.values
        new_spreadsheet_df = cosine_matrix_df.loc[user_gene_list].sum()
        property_rank_df[col_name] = new_spreadsheet_df.values

    cosine_sum_val = np.ravel(property_rank_df.values).round(12)
    set_name = np.array(list(property_rank_df.columns.values) * (property_rank_df.shape[0]))
    gene_name = np.repeat(property_rank_df.index.values, property_rank_df.shape[1])

    ret_col = ['user_gene_set', 'property_gene_set', 'cosine_sum']
    result_val = np.column_stack((set_name, gene_name, cosine_sum_val))
    result_df = pd.DataFrame(result_val, columns=ret_col).sort_values("cosine_sum", ascending=0)
    return result_df


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


# A global list used exclusively in get_fisher_exact_test as a callback variable of parallel execution
fisher_contingency_pval_parallel_insertion = []


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

    dimension = [range(overlap_count.shape[0]), range(overlap_count.shape[1])]
    combinations = list(itertools.product(*dimension))
    parallelism = dstutil.determine_parallelism_locally(len(combinations))

    try:
        p = multiprocessing.Pool(processes=parallelism)
        p.starmap_async(fisher_exact_worker, zip(itertools.repeat(sparse_dict),
                                                 itertools.repeat(overlap_count),
                                                 itertools.repeat(user_count),
                                                 itertools.repeat(gene_count),
                                                 itertools.repeat(universe_count),
                                                 itertools.repeat(set_list),
                                                 combinations),
                        callback=callback_extend_list)
        p.close()
        p.join()
        # print(fisher_contingency_pval_parallel_insertion)
        # print(type(fisher_contingency_pval_parallel_insertion))
        return fisher_contingency_pval_parallel_insertion
    except:
        raise OSError("Failed running parallel processing:{}".format(sys.exc_info()))


def callback_extend_list(item):
    """ Callback function called by function get_fisher_exact_test. It extends items to a global list

    Args:
        item: input to append to global list

    Returns:
         None
    """
    fisher_contingency_pval_parallel_insertion.extend(item)


def fisher_exact_worker(sparse_dict, overlap_count, user_count, gene_count, universe_count, set_list, combinations):
    """ worker of parallel execution called by function get_fisher_exact_test

    Args:
        sparse_dict:
        overlap_count:
        user_count:
        gene_count:
        universe_count:
        set_list:
        combinations:

    Returns:
        row_item
    """
    i, j = combinations[0], combinations[1]
    table = build_fisher_contingency_table(overlap_count[i, j], user_count[j], gene_count[0, i], universe_count)
    pvalue = stats.fisher_exact(table, alternative="greater")[1]
    new_pval = np.round(-1.0 * np.log10(pvalue), 12)
    row_item = [set_list[j], sparse_dict[i], new_pval, int(universe_count), 
    int(user_count[j]), int(gene_count[0, i]), int(overlap_count[i, j])]
    return row_item


def save_fisher_test_result(fisher_contingency_pval, results_dir, set_list, threshold):
    """ Save two output files of fisher exact test results.

    Args:
        fisher_contingency_pval: list of seven items lists.
        set_list: column values of spreadsheet.
    Returns:
        result_df: the final dataframe of fisher exact test
        threshold: only return the pvalues above the threshold
    """
    df_col = ["user_gene_set", "property_gene_set", "pval", "universe_count", \
              "user_count", "property_count", "overlap_count"]
    result_df = pd.DataFrame(
        fisher_contingency_pval, columns=df_col).sort_values("pval", ascending=0)

    result_df_with_score = pd.DataFrame(columns=set_list)
    for gene_set in set_list:
        result_df_with_score.loc[:, gene_set] = result_df[result_df['user_gene_set'] == gene_set].values[:, 1]
    save_timestamped_df(result_df_with_score, results_dir, 'fisher_ranked_by_property')

    result_df_with_rank = result_df
    if len(fisher_contingency_pval) > 100:
        result_df_with_rank = result_df[result_df['pval'] > threshold]
    
    save_timestamped_df(result_df_with_rank, results_dir, 'fisher_sorted_by_property_score')

    return result_df


def rank_drawr_property(final_spreadsheet_df, pg_network_n1_names):
    """ This is to rank properties for each user gene set.

    Args:
        final_spreadsheet_df: final smoothed dataframe.
        pg_network_n1_names: list of property names.
    Returns:
        prop_spreadsheet_df: a new spreadsheet with user gene set as header and properties as values.
    """
    prop_spreadsheet_df = final_spreadsheet_df.loc[pg_network_n1_names]
    prop_spreadsheet_df.iloc[:, :-1] = prop_spreadsheet_df.iloc[:, :-1].apply(
        lambda x: (x - prop_spreadsheet_df['base']).sort_values(ascending=0).index.values)
    prop_spreadsheet_df = prop_spreadsheet_df.drop('base', 1)

    return prop_spreadsheet_df


def construct_drawr_result_df(input_df, start_index, end_index, map_back, run_parameters):
    """Construct a five-column DRaWR result dataframe with
    selected rows from smoothed spreadsheet dataframe.
    Args:
        input_df: input spreadsheet dataframe.
        start_index: starting index
        end_index: end index
        map_back: boolean value to check whether to map gene names to the orignal names
    Returns:
        result_df: result five-column dataframe.
    """
    len_set_names = input_df.shape[1] - 1
    smooth_base = input_df.values[start_index:end_index, -1]
    smooth_base = smooth_base[:, np.newaxis]
    diff_smooth_spreadsheet = input_df.values[start_index:end_index, :-1] - smooth_base
    diff_smooth_spreadsheet /= np.abs(np.max(diff_smooth_spreadsheet, axis=0))

    diff_val = np.ravel(diff_smooth_spreadsheet)
    orig_val = np.ravel(input_df.values[start_index:end_index, :-1])
    set_name = np.array(list(input_df.columns.values[:-1]) * (end_index - start_index))
    input_gene_name = input_df.index.values[start_index:end_index]

    if map_back is True:
        map_df = pd.read_csv(run_parameters["gene_names_map"], index_col=0, header=None, sep='\t')
        input_gene_name = map_df.loc[input_gene_name].values
        ret_col = ['user_gene_set', 'gene_node_id', 'difference_score', 'query_score', 'baseline_score']
    else:
        ret_col = ['user_gene_set', 'property_gene_set', 'difference_score', 'query_score', 'baseline_score']
    new_gene_name = np.repeat(input_gene_name, len_set_names)
    base_val = np.repeat(input_df['base'].values[start_index:end_index], len_set_names)
    
    result_val = np.column_stack((set_name, new_gene_name, diff_val, orig_val, base_val))
    result_df = pd.DataFrame(result_val, columns=ret_col).sort_values("difference_score", ascending=0)
    result_df = result_df[result_df['difference_score'] > 0.5]
    return result_df


def save_timestamped_df(input_df, results_dir, output_file_name):
    """ Save dataframe to files with timestamped name.

    Args:
        fisher_contingency_pval: list of seven items lists.
        results_dir: directory to save outputs.
        output_file_name: file name.
    """
    file_name = kn.create_timestamped_filename(output_file_name, "df")
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
    map_df = pd.read_csv(run_parameters["gene_names_map"], index_col=0, header=None, sep='\t')
    new_droplist_df = pd.DataFrame(map_df.loc[droplist].values, columns=[droplist_name])
    file_name = kn.create_timestamped_filename(droplist_name, "tsv")
    kn.save_df(new_droplist_df, run_parameters['results_directory'], file_name)


def build_hybrid_sparse_matrix(run_parameters, normalize_by_sum, construct_by_union):
    """This is to build hybrid sparse matrix with gene gene network and
    gene property network.

    Args:
        run_parameters: dictionary of run parameters.
        normalize_by_sum: boolean value to check normalization.
        construct_by_union: boolean value to check construct by union.

    Returns:
        network_sparse: output sparse matrix.
        unique_gene_names: gene names of the hybrid matrix
        pg_network_n1_names: property names of hybrid matrix.
    """
    pg_network_df = kn.get_network_df(run_parameters['pg_network_name_full_path'])
    gg_network_df = kn.get_network_df(run_parameters['gg_network_name_full_path'])

    pg_network_n1_names, \
    pg_network_n2_names = kn.extract_network_node_names(pg_network_df)

    gg_network_n1_names, \
    gg_network_n2_names = kn.extract_network_node_names(gg_network_df)

    # limit the gene set to the intersection of networks (gene_gene and prop_gene) and user gene set
    unique_gene_names = kn.find_unique_node_names(gg_network_n1_names, gg_network_n2_names)

    if construct_by_union is True:
        unique_gene_names = kn.find_unique_node_names(unique_gene_names, pg_network_n2_names)
    else:
        pg_network_df = kn.update_network_df(pg_network_df, unique_gene_names, 'node_2')

    unique_gene_names_dict = kn.create_node_names_dict(unique_gene_names)
    pg_network_n1_names_dict = kn.create_node_names_dict(
        pg_network_n1_names, len(unique_gene_names))

    unique_all_node_names = unique_gene_names + pg_network_n1_names
    # map every gene name to a sequential integer index
    gg_network_df = kn.map_node_names_to_index(gg_network_df, unique_gene_names_dict, "node_1")
    gg_network_df = kn.map_node_names_to_index(gg_network_df, unique_gene_names_dict, "node_2")
    pg_network_df = kn.map_node_names_to_index(pg_network_df, pg_network_n1_names_dict, "node_1")
    pg_network_df = kn.map_node_names_to_index(pg_network_df, unique_gene_names_dict, "node_2")

    gg_network_df = kn.symmetrize_df(gg_network_df)
    pg_network_df = kn.symmetrize_df(pg_network_df)

    if normalize_by_sum is True:
        gg_network_df = kn.normalize_network_df_by_sum(gg_network_df, 'wt')
        pg_network_df = kn.normalize_network_df_by_sum(pg_network_df, 'wt')

    hybrid_network_df = kn.form_hybrid_network_df([gg_network_df, pg_network_df])

    # store the network in a csr sparse format
    network_sparse = kn.convert_network_df_to_sparse(
        hybrid_network_df, len(unique_all_node_names), len(unique_all_node_names))

    return network_sparse, unique_gene_names, pg_network_n1_names
