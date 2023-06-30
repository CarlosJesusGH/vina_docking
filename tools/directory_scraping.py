import glob, os

def list_files_in_dir(directory, search_pattern='*'):
    # list files in directories
    file_paths = sorted(glob.glob(directory + search_pattern))
    file_names = list(map(lambda x: os.path.splitext(os.path.basename(x))[0], file_paths))
    return file_paths, file_names

def list_files_in_dir_recursively(directory, search_pattern='*'):
    # list files in directories
    file_paths = sorted(glob.glob(directory + "/**/" + search_pattern, recursive=True))
    # file_names = list(map(lambda x: os.path.splitext(os.path.basename(x))[0], file_paths))
    file_names = list(map(lambda x: os.path.basename(x), file_paths))
    subdirectory_names = list(map(lambda x: os.path.basename(os.path.dirname(x)), file_paths))
    return file_paths, file_names, subdirectory_names

def mkdir_if_not_exists(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

def mkdir_or_clear(directory):
    rmdir_if_exists(directory)
    mkdir_if_not_exists(directory)

def rmdir_if_exists(directory):
    # import shutil
    if os.path.exists(directory):
        # shutil.rmtree(directory)
        os.system("rm -rf " + directory)
        
def prepare_directory_from_dict(dirs_dict):
    for d in dirs_dict:
        if not os.path.exists(dirs_dict[d]):
            os.makedirs(dirs_dict[d])