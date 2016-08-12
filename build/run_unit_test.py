import subprocess

def run_unit_test():

    cmd_pre = ['python3', '../test/unit/test_perform_fisher_exact_test.py']
    p_pre = subprocess.Popen(cmd_pre, stdout=subprocess.PIPE, shell=False)
    p_pre.wait()


run_unit_test()

