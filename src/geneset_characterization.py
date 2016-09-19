"""
Created on Wed Jul 20 14:47:45 2016
@author: The Gene Sets Characterization dev team
"""
import sys
import time
import yaml
import os   

def fisher(run_parameters):
    """fisher geneset characterization

    Args:
        run_parameters: parameter set dictionary.
    """
    import geneset_characterization_toolbox as tl
    t0 = time.time()
    tl.run_fisher(run_parameters)
    return time.time() - t0
    
def DRaWR(run_parameters):
    """Discriminative Random Walk with Restart

    Args:
        run_parameters: parameter set dictionary.
    """
    import geneset_characterization_toolbox as tl
    t0 = time.time()
    tl.run_DRaWR(run_parameters)
    return time.time() - t0

def net_one(run_parameters):
    """net_one geneset characterization method

    Args:
        run_parameters: parameter set dictionary.
    """
    import geneset_characterization_toolbox as tl
    t0 = time.time()
    tl.run_net_one(run_parameters)
    return time.time() - t0

SELECT = {
    "fisher":fisher,
    "DRaWR":DRaWR,
    "net_one":net_one}

def main():
    """
    This the main function to characterize gene
    set from three methods: fisher, DRaWR, net one.
    """

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
