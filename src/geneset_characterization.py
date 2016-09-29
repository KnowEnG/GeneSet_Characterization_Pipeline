"""
Created on Wed Jul 20 14:47:45 2016
@author: The Gene Sets Characterization dev team
"""

def fisher(run_parameters):
    """fisher geneset characterization

    Args:
        run_parameters: parameter set dictionary.
    """
    import geneset_characterization_toolbox as tl
    tl.run_fisher(run_parameters)
    
def DRaWR(run_parameters):
    """Discriminative Random Walk with Restart

    Args:
        run_parameters: parameter set dictionary.
    """
    import geneset_characterization_toolbox as tl
    tl.run_DRaWR(run_parameters)

def net_one(run_parameters):
    """net_one geneset characterization method

    Args:
        run_parameters: parameter set dictionary.
    """
    import geneset_characterization_toolbox as tl
    tl.run_net_one(run_parameters)

SELECT = {
    "fisher":fisher,
    "DRaWR":DRaWR,
    "net_one":net_one}

def main():
    """
    This the main function to characterize gene
    set from three methods: fisher, DRaWR, net one.
    """
    import yaml
    from knpackage.toolbox import get_run_directory_and_file
    from knpackage.toolbox import get_run_parameters

    run_directory, run_file = get_run_directory_and_file(sys.argv)
    run_parameters = get_run_parameters(run_directory, run_file)
    SELECT[run_parameters["method"]](run_parameters)
    

if __name__ == "__main__":
    main()
