"""
Created on Wed Jul 20 14:47:45 2016
@author: The Gene Sets Characterization dev team
"""

def fisher(run_parameters):
    '''fisher geneset characterization'''
    from geneset_characterization_toolbox import run_fisher
    run_fisher(run_parameters)
    
def DRaWR(run_parameters):
    '''Discriminative Random Walk with Restart'''
    from geneset_characterization_toolbox import run_DRaWR
    run_DRaWR(run_parameters)

def net_path(run_parameters):
    '''net_path geneset characterization method'''
    from geneset_characterization_toolbox import run_net_path
    run_net_path(run_parameters)

SELECT = {
    "fisher":fisher,
    "DRaWR":DRaWR,
    "net_path":net_path}

def main():
    """
    This the main function to characterize gene
    set from three methods: fisher, DRaWR, net path.
    """
    import sys
    from knpackage.toolbox import get_run_directory_and_file
    from knpackage.toolbox import get_run_parameters

    run_directory, run_file = get_run_directory_and_file(sys.argv)
    run_parameters = get_run_parameters(run_directory, run_file)
    SELECT[run_parameters["method"]](run_parameters)
    

if __name__ == "__main__":
    main()
