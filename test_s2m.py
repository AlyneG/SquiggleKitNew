import os
import sys
from ont_fast5_api.conversion_tools import single_to_multi_fast5
from ont_fast5_api.conversion_tools.single_to_multi_fast5 import batch_convert_single_to_multi
from ont_fast5_api.conversion_tools.single_to_multi_fast5 import create_multi_read_file
# batch_conver_single_to_multi(input_path, output_folder, filename_base, batch_size, threads, recursive, follow_symlinks, target_compression)
from ont_fast5_api.multi_fast5 import MultiFast5File
from ont_fast5_api.multi_fast5 import Fast5File


def main():
    s2m('./fast5', './fast5/output.fast5', None)


def s2m(f5_path, output_file, target_compression):
    '''
    Combine single fast5 files from one directory into 1 multi fast5 file
    '''
    # batch_convert_single_to_multi(input_path, output_folder, filename_base, batch_size, threads, recursive, follow_symlinks, target_compression)
    # create_multi_read_file(input_files, output_file, target_compression)
    # target_compression (str) [vbz, gzip], target compresesion format

    filenames = []
    results = []
    output = None
    for dirpath, dirnames, files in os.walk(f5_path):
        for fast5 in files:
            if fast5.endswith('.fast5'):
                filenames.append(os.path.join(dirpath, fast5))
    if filenames:
        results, output = create_multi_read_file(filenames, output_file, target_compression)
    
    print(filenames, results, output)


if __name__ == '__main__':
    main()
