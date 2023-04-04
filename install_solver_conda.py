import tarfile
import os
import sys
import shutil

root_path = os.path.abspath(os.path.dirname(__file__))
lib_a_map = {'minos56': 'minos', 'qminos56': 'quadminos'}
def install_qminos(solver_path: str = root_path + 'solver/'):
    """Create the Python qminos solver in the environment. Run with Makefile. 
    Qminos must be downloaded

    Parameters
    ----------
    solver_path : str
        by default, assumes the solver directory has been made in the package
    """
    # STEP 1: get the qminos solver locally
    # TODO: replace STEP1 with juts downloading qminos if possible -- can actually combine with human_me/data/_download_data.py (solver and data subdirs all into same parent directory)
    if not os.path.isdir(solver_path):
        raise ValueError('The path to the speciifed olver directory does not exist')
    if not os.path.isdir(os.path.join(solver_path, 'qminos1114')):
        if not os.path.isfile(os.path.join(solver_path, 'qminos1114b.tar.gz')):
            raise ValueError('The qminos solver must be downloaded and named "qminos1114" unzipped or "qminos1114b.tar.gz" zipped')
        else:
            with tarfile.open(os.path.join(solver_path, 'qminos1114b.tar.gz'), mode='r|gz') as tar_f:
                tar_f.extractall(path=solver_path)
    # STEP2 : get the solvemepy github repo locally
    os.system(' '.join(['git -C', solver_path, 'clone https://github.com/SBRG/solvemepy.git'])) #@2a2c9c098d5bad957ef41637955fe338a31bac4c
    # STEP 3: copy makefiles within qminos and make 
    # steps 3-4 as explained in https://github.com/SBRG/solvemepy
    for dir_m in ['minos56', 'qminos56']:
        shutil.copy2(src=os.path.join(solver_path, 'qminos1114', 'Makefile.defs'),
                    dst=os.path.join(solver_path, 'qminos1114', dir_m))
        os.system(' '.join(['make clean -C', os.path.join(solver_path, 'qminos1114', dir_m)]))
        os.system(' '.join(['make -C', os.path.join(solver_path, 'qminos1114', dir_m)]))

        # STEP 4: transfer files from qminos to solvemepy 
        shutil.copy2(src=os.path.join(solver_path, 'qminos1114', dir_m, 'lib/lib' + lib_a_map[dir_m] + '.a'), 
                    dst=os.path.join(solver_path, 'solvemepy'))

    os.chdir(os.path.join(solver_path, 'solvemepy'))
    os.system('python ' + os.path.join(solver_path, 'solvemepy/setup.py') + ' develop')
    os.chdir(root_path)

def install_solver_(solver_path: str = root_path + 'solver/', solver_type = 'qminos'):
    if solver_type == 'qminos':
        install_qminos(solver_path)
    else:
        raise ValueError('Only the qminos solver can be installed for now')

if __name__ == "__main__":
    install_solver_(sys.argv[1], sys.argv[2])
