

def batch_download_receptors(list_receptors_file, raw_receptors_dir):
    import subprocess as sp
    output = sp.getoutput('/home/jovyan/vina_docking/tools/3dparty_scripts/batch_download.sh -f ' + list_receptors_file + ' -o ' + raw_receptors_dir + ' -p ')
    return output