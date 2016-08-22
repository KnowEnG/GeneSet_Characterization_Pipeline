"""
Created on Wed Jul 20 14:47:45 2016
@author: Xi, Suyang, Jing, Dan, Nahil
"""

import time
import yaml
import os    
    
def fisher(run_parameters):
    '''fisher geneset characterization'''
    import geneset_characterization_toolbox as tl
    t0 = time.time()
    tl.run_fisher(run_parameters)
    # print('run_fisher total time: {}'.format(time.time() - t0))
    return time.time() - t0

def DRaWR(run_parameters):
    '''Discriminative Random Walk with Restart'''
    import geneset_characterization_toolbox as tl
    t0 = time.time()
    tl.run_DRaWR(run_parameters)
    # print('run_DRaWR total time: {}'.format(time.time() - t0))
    return time.time() - t0

def net_one(run_parameters):
    '''net_one geneset characterization method'''
    import geneset_characterization_toolbox as tl
    t0 = time.time()
    tl.run_net_one(run_parameters)
    # print('run_net_one total time: {}'.format(time.time() - t0))
    return time.time() - t0

SELECT = {
    "fisher":fisher,
    "DRaWR":DRaWR,
    "net_one":net_one}

def main():
    import sys
    from knpackage.toolbox import get_run_directory_and_file
    from knpackage.toolbox import get_run_parameters

    run_directory, run_file = get_run_directory_and_file(sys.argv)
    run_parameters = get_run_parameters(run_directory, run_file)

    f_time = SELECT[run_parameters["method"]](run_parameters)
    run_parameters["running_time"] = f_time
    
    with open(os.path.join(run_directory, run_file), 'w') as file_handle:
        yaml.dump(run_parameters, file_handle, default_flow_style=False)

if __name__ == "__main__":
    main()