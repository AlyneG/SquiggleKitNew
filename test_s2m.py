import os
import sys
from ont_fast5_api.conversion_tools import single_to_multi_fast5
from ont_fast5_api.conversion_tools.single_to_multi_fast5 import create_multi_read_file
# batch_conver_single_to_multi(input_path, output_folder, filename_base, batch_size, threads, recursive, follow_symlinks, target_compression)
from ont_fast5_api.multi_fast5 import MultiFast5File
from ont_fast5_api.multi_fast5 import Fast5File
import shutil


def main():
    os.mkdir("fast5_fetcher_temp")
    tmp_path = os.path.abspath("fast5_fetcher_temp")
    s2m('./fast5', tmp_path, 'output.fast5', None)
    s2m('./fast5', tmp_path, 'output2.fast5', None)

def s2m(f5_path, save_path, output_file, target_compression):
    '''
    Combine single fast5 files into 1 multi fast5 file
    '''
    filenames = []
    results = []
    output = None
    
    for dirpath, dirnames, files in os.walk(f5_path):
        for fast5 in files:
            if fast5.endswith('.fast5'):
                filenames.append(os.path.join(dirpath, fast5))
    if filenames:
        results, output = create_multi_read_file(filenames, os.path.join(save_path, output_file), target_compression)
    
    with open(os.path.join(save_path, "filename_mapping.txt"), 'a') as out_sum:
        # get readIDs from sequencing summary?
        for S in results:
           out_sum.write("{}\t{}\n".format(output_file, S))


if __name__ == '__main__':
    main()
