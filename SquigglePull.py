import os
import sys
import argparse
import traceback
import numpy as np
import h5py
import time
import sklearn.preprocessing

'''
    SquigglePull
    James M. Ferguson (j.ferguson@garvan.org.au)
    Genomic Technologies
    Garvan Institute
    Copyright 2017

    Pull squiggle data from fast5 files

    input:
        - path to fast5 files

    output:
        - tsv signal file

    TODO:
        - Dynamic columns and data types
        - Multi fast5 file support
        - paf, sam, fastq, or flat file support for filtering
        - multiprocessing
        - use # notation at start of file for static values, size reduction,
          including things like kits, flowcells, versions, etc, for comparisons.



    Testing:
        python SquigglePull.py -es -p test/R9_event_data/ -t 20,110 -f pos1 > data.tsv
        python SquigglePull.py -r -p test/R9_raw_data/ -f all > data.tsv

    Notes:
        should do some target-type validations before executing and exit.
    -----------------------------------------------------------------------------
    MIT License

    Copyright (c) 2017 James Ferguson

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.

'''


class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)


def main():
    '''
    One function to rule them all, one function to find them, One function to bring them
    all and in the darkness call them.
    '''

    parser = MyParser(
        description="SquigglePull - extraction and (optional) conversion to pA of raw signal from Oxford Nanopore fast5 files")

    # arguments
    parser.add_argument("-p", "--path",
                        help="Top directory path of fast5 files")
    parser.add_argument("--single", action="store_true",
                        help="single_fast5 files")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Engage higher output verbosity")
    parser.add_argument("-s", "--scale", action="store_true",
                        help="Scale signal output for comparison")
    parser.add_argument("-r", "--raw_signal", action="store_true",
                        help="No conversion to pA, raw signal is extracted instead")
    #parser.add_argument("--round", help="Value for num of decimal places converted signal is to be rounded to")
    parser.add_argument("-i", "--extra_info", action="store_true",
                        help="Print extra information used for signal conversion and in methylation calling - nanopolish/f5c")
    # parser.add_argument("-a", "--paf",
    #                     help="paf alignment file for nt approach - Benchmarking")
    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    if args.verbose:
        sys.stderr.write("Verbose mode on. Starting timer")
        start_time = time.time()


    # process fast5 files given top level path
    # This should work for multi-fast5 too, push detect into extract_f5()
    for dirpath, dirnames, files in os.walk(args.path):
        for fast5 in files:
            if fast5.endswith('.fast5'):
                fast5_file = os.path.join(dirpath, fast5)
                # extract data from file
                data = extract_f5_all(fast5_file, args)
                if not data:
                    sys.stderr.write("main():data not extracted from {}. Moving to next file.".format(fast5_file))
                    continue
                # print data
                if args.single:
                    print_data(data, args, fast5)
                else:
                    for read in data:
                        print_data(data[read], args, fast5)
    if args.verbose:
        end_time = time.time() - start_time
        sys.stderr.write("Time taken: {}\n".format(end_time))


#def extract_f5(filename, args, batch=False):
    '''
    inputs:
        filepath/name
        optional:
            Raw vs Events
            batch flags
    does:
        open fast5 files, extract whole signal and read data
    Returns:
        dic for further processing

    2 methods:
        - Open one at a time (slow) - single thread
        - Open batches at a time (fast) - paralellised


    Takes the latest basecalled events table.
    '''

#    f5_dic = {'raw': [], 'seq': '', 'readID': '',
#            'digitisation': 0.0, 'offset': 0.0, 'range': 0.0,
#            'sampling_rate': 0.0}

    # open fast5 file
#    try:
#        hdf = h5py.File(filename, 'r')
#    except:
#        traceback.print_exc()
#        sys.stderr.write("extract_fast5():fast5 file failed to open: {}".format(filename))
#        f5_dic = {}
#        return f5_dic

    # extract information
#    try:
#        c = list(hdf['Raw/Reads'].keys())
#        for col in hdf['Raw/Reads/'][c[0]]['Signal'][()]:
#            f5_dic['raw'].append(int(col))
#        f5_dic['readID'] = hdf['Raw/Reads/'][c[0]].attrs['read_id'].decode()
#        if not(args.raw_signal) or args.extra_info:
#            f5_dic['digitisation'] = hdf['UniqueGlobalKey/channel_id'].attrs['digitisation']
#            f5_dic['offset'] = hdf['UniqueGlobalKey/channel_id'].attrs['offset']
#            f5_dic['range'] = float("{0:sud.2f}".format(hdf['UniqueGlobalKey/channel_id'].attrs['range']))
#            f5_dic['sampling_rate'] = hdf['UniqueGlobalKey/channel_id'].attrs['sampling_rate']
#    except:
#        traceback.print_exc()
#        sys.stderr.write("extract_fast5():failed to extract raw signal or fastq from {}".format(filename))
#        f5_dic = {}

#    return f5_dic

#def read_multi_fast5(filename, args):
    '''
    read multifast5 file and return data
    '''
#    f5_dic = {}
#    with h5py.File(filename, 'r') as hdf:
#        for read in list(hdf.keys()):
#            f5_dic[read] = {'raw': [], 'seq': '', 'readID': '', 'digitisation': 0.0,
#                            'offset': 0.0, 'range': 0.0, 'sampling_rate': 0.0}

#            try:
#                for col in hdf[read]['Raw/Signal'][()]:
#                    f5_dic[read]['raw'].append(int(col))
#                f5_dic[read]['readID'] = hdf[read]['Raw'].attrs['read_id'].decode()

#                if not(args.raw_signal) or args.extra_info:
#                    f5_dic[read]['digitisation'] = hdf[read]['channel_id'].attrs['digitisation']
#                    f5_dic[read]['offset'] = hdf[read]['channel_id'].attrs['offset']
#                    f5_dic[read]['range'] = float("{0:.2f}".format(hdf[read]['channel_id'].attrs['range']))
#                    f5_dic[read]['sampling_rate'] = hdf[read]['channel_id'].attrs['sampling_rate']
#                
#            except:
#                traceback.print_exc()
#                sys.stderr.write("extract_fast5():failed to read readID: {}".format(read))
#    return f5_dic

def extract_f5_all(filename, args):
    '''
    inputs:
        filepath/name
        args from command line
    does:
        open fast5 files, extract whole signal and read data and converts to pA by default
    Returns:
        dic for further processing/printing
    '''
    f5_dic = {}
    with h5py.File(filename, 'r') as hdf:
        # single fast5 files
        if args.single:
            f5_dic = {'raw': [], 'seq': '', 'readID': '',
                    'digitisation': 0.0, 'offset': 0.0, 'range': 0.0,
                    'sampling_rate': 0.0}
            # extract the data
            try:
                c = list(hdf['Raw/Reads'].keys())
                for col in hdf['Raw/Reads/'][c[0]]['Signal'][()]:
                    f5_dic['raw'].append(int(col))
                f5_dic['readID'] = hdf['Raw/Reads/'][c[0]].attrs['read_id'].decode()
                digitisation = hdf['UniqueGlobalKey/channel_id'].attrs['digitisation']
                offset = hdf['UniqueGlobalKey/channel_id'].attrs['offset']
                range = float("{0:.2f}".format(hdf['UniqueGlobalKey/channel_id'].attrs['range']))
                
                # convert to pA                    
                if not(args.raw_signal):
                    f5_dic['raw'] = np.array(f5_dic['raw'], dtype=int)
                    f5_dic['raw'] = convert_to_pA_numpy(f5_dic['raw'], digitisation, range, offset)
                    f5_dic['raw'] = np.round(f5_dic['raw'], 2)

                # save the extra info for printing                
                if args.extra_info:
                    f5_dic['digitisation'] = digitisation
                    f5_dic['offset'] = offset
                    f5_dic['range'] = range
                    f5_dic['sampling_rate'] = hdf['UniqueGlobalKey/channel_id'].attrs['sampling_rate']
            except:
                traceback.print_exc()
                sys.stderr.write("extract_fast5_all():failed to extract raw signal or fastq from {}".format(filename))
                f5_dic = {}

        # multi fast5 files
        else:
            for read in list(hdf.keys()):
                f5_dic[read] = {'raw': [], 'seq': '', 'readID': '', 
                                'digitisation': 0.0, 'offset': 0.0, 'range': 0.0,
                                'sampling_rate': 0.0}

                # extract the data
                try:
                    for col in hdf[read]['Raw/Signal'][()]:
                        f5_dic[read]['raw'].append(int(col))
                    f5_dic[read]['readID'] = hdf[read]['Raw'].attrs['read_id'].decode()
                    digitisation = hdf[read]['channel_id'].attrs['digitisation']
                    offset = hdf[read]['channel_id'].attrs['offset']
                    range = float("{0:.2f}".format(hdf[read]['channel_id'].attrs['range']))

                    # convert to pA
                    if not(args.raw_signal):
                        f5_dic[read]['raw'] = np.array(f5_dic[read]['raw'], dtype=int)
                        f5_dic[read]['raw'] = convert_to_pA_numpy(f5_dic[read]['raw'], digitisation, range, offset)
                        f5_dic[read]['raw'] = np.round(f5_dic[read]['raw'], 2)
                    
                    # save the extra info for printing                    
                    if args.extra_info:
                        f5_dic[read]['digitisation'] = digitisation
                        f5_dic[read]['offset'] = offset
                        f5_dic[read]['range'] = range
                        f5_dic[read]['sampling_rate'] = hdf[read]['channel_id'].attrs['sampling_rate']
                    
                except:
                    traceback.print_exc()
                    sys.stderr.write("extract_fast5_all():failed to read readID: {}".format(read))
    return f5_dic

'''
def pull_target(data, args, min_length=50, paf=None):
    
    Pull out target region from data.

    inputs:
        - data - dictionary containing reads
        - target - pos1: 20,110 - event/raw positions
        - target_type - pos1

    does:
        ...explain methods...

    Returns:
        - Regions of interest labelled by read_id/filename

    dicf5_dic = {'events': [], 'moves': [], 'seq': '', 'readID': ''}
    
    default = []
    region = []

    signal = np.array(data['raw'])
    if args.scale:
        signal = scale_data(signal)
    target = str(len(signal))

    region.append(target)
    region.append(signal)

    if region:
        return region
    else:
        sys.stderr.write("pull_target():Something went wrong. Region not found")
        return default
'''

def scale_data(data):
    '''
    Scale shift and scale for comparisons
    '''
    try:
        scaled_data = sklearn.preprocessing.scale(data,
                                                  axis=0,
                                                  with_mean=True,
                                                  with_std=True,
                                                  copy=True)
    except:
        traceback.print_exc()
        sys.stderr.write("scale_data():Something went wrong, failed to scale data")
        return 0
    return scaled_data

#def convert_to_pA(d):
    '''
    convert raw signal data to pA using digitisation, offset, and range
    float raw_unit = range / digitisation;
    for (int32_t j = 0; j < nsample; j++) {
        rawptr[j] = (rawptr[j] + offset) * raw_unit;
    }
    '''
#    digitisation = d['digitisation']
#    range = d['range']
#    offset = d['offset']
#    raw_unit = range / digitisation
#    new_raw = []
#    for i in d['raw']:
#        j = (i + offset) * raw_unit
#        new_raw.append("{0:.2f}".format(round(j,2)))
#    return new_raw

def convert_to_pA_numpy(d, digitisation, range, offset):
    raw_unit = range / digitisation
    return (d + offset) * raw_unit

def print_data(data, args, fast5):
    if args.scale:
        data['raw'] = scale_data(data['raw'])

    ar = []
    for i in data['raw']:
        ar.append(str(i))

    if args.extra_info:
        print('{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(fast5, data['readID'],
                data['digitisation'], data['offset'], data['range'],
                data['sampling_rate'], '\t'.join(ar)))
    else:                        
        print('{}\t{}\t{}'.format(
                fast5, data['readID'], '\t'.join(ar)))

if __name__ == '__main__':
    main()

