"""
lanier4@illinois.edu

"""
import os
import filecmp
import time

GSC_options_dict = {    'run_fisher'        : 'BENCHMARK_1_fisher',
                        'run_drawr'         : 'BENCHMARK_2_DRaWR',
                        'run_netpath_small' : 'TEST_3_net_path'     }

verif_dir = '../data/verification'
results_dir = './run_dir/results'

file_ext = '.tsv'

def run_all_BENCHMARKs_and_TESTs():
    """ run the make file targes for all yaml files and compre the results with their verification files """
    print('\n\nBegin Verification Testing:\t\t', time.strftime("%a, %b %d, %Y at %H:%M:%S", time.localtime()))
    directory_methods_dict = {v: k for k, v in (GSC_options_dict).items()}
    verification_directory_list = sorted(directory_methods_dict.keys())

    for test_directory in verification_directory_list:
        verification_directory = os.path.join(verif_dir, test_directory)
        verification_method = 'make' + ' ' + directory_methods_dict[test_directory]
        print('\n\n\n\tRun Method:', verification_method, '\n', verification_directory)

        os.system(verification_method)
        mismatch_list = python_dir_compare(verification_directory, results_dir)
        if len(mismatch_list) > 0:
            print('\t\t\t\t Files Not Verified \t\t FAIL')
            for missed_file in mismatch_list:
                print('\n\t\t FILES NOT VERIFIED:\n', missed_file, '\n')
        else:
            print('\t\t\t\t All Files Verified \t\t PASS')

        print('removing result files:')
        for tmp_file_name in os.listdir(results_dir):
            if os.path.isfile(os.path.join(results_dir, tmp_file_name)):
                os.remove(os.path.join(results_dir, tmp_file_name))

    print('\n\nFinished Verification Testing:\t\t', time.strftime("%a, %b %d, %Y at %H:%M:%S", time.localtime()))

def python_dir_compare(verif_dir, results_dir):
    mismatch_list = []
    veri_files_list = os.listdir(verif_dir)
    res_files_list = os.listdir(results_dir)

    for veri_file in veri_files_list:
        if os.path.isfile(os.path.join(verif_dir, veri_file)) == True:
            if veri_file[-len(file_ext):] == file_ext:
                file_pre = veri_file[0:-len(file_ext)]
                for res_file in res_files_list:
                    if os.path.isfile(os.path.join(results_dir, res_file)) and res_file[0:len(file_pre)] == file_pre:
                        res_file = os.path.join(results_dir, res_file)
                        FILES_MATCH = filecmp.cmp(os.path.join(verif_dir, veri_file), res_file)
                        if FILES_MATCH == False:
                            mismatch_list.append(res_file)
    return mismatch_list

def main():
    try:
        print('environment setup:')
        os.system('make env_setup')
    except:
        print('environment setup: FAIL')
        pass

    run_all_BENCHMARKs_and_TESTs()
    print('\n')

if __name__ == "__main__":
    main()