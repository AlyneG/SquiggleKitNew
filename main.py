from flask import Flask, render_template, request, redirect, Markup
import os
import sys
import time
import h5py
import traceback
import numpy as np
#from bokeh.plotting import figure
#from bokeh.client import pull_session
#from bokeh.embed import server_session
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 14})
import mpld3
from mpld3 import plugins
import json

app = Flask(__name__)

@app.route("/",methods = ["POST", "GET"])
def home():
    if request.method == "POST":
        f5_path = request.form.get('f5_path')
        type = request.form.get('type')
        if not os.path.isdir(f5_path):
            return render_template("home.html", error=True, f5_path=f5_path)
        else:
            return render_template("loading.html", f5_path=f5_path, type=type)
            
    return render_template("home.html")

@app.route("/results")
def results(f5_path=None):

    f5_path = request.args['f5_path']
    type = request.args['type']

    if request.args['processing'] == '0':
        if os.path.isfile(f5_path+"/data_"+type+".tsv"):
            exception = "Data file for already exists"
            return render_template("exception.html", f5_path=f5_path, type=type, exception=exception)

    if request.args['processing'] == '1':
        os.remove(f5_path+"/data_"+type+".tsv")

    count = 0
    with open(os.path.join(f5_path, "data_"+type+".tsv"), 'a') as out_sum:
        for dirpath, dirnames, files in os.walk(f5_path):
            for fast5 in files:
                if fast5.endswith('.fast5'):
                    fast5_file = os.path.join(dirpath, fast5)
                    # extract data from file
                    data, multi = extract_f5_all(fast5_file, request.args['type'])
                    #print data to a single file
                    if not multi:
                        count += 1
                        ar = []
                        for i in data['raw']:
                            ar.append(str(i))

                        out_sum.write('{}\t{}\t{}\n'.format(
                                fast5, data['readID'], '\t'.join(ar)))
                    else:
                        for read in data:
                            count += 1
                            ar = []
                            for i in data[read]['raw']:
                                ar.append(str(i))

                            out_sum.write('{}\t{}\t{}\n'.format(
                                    fast5, data[read]['readID'], '\t'.join(ar)))
    return render_template("results.html", f5_path=f5_path, type=type, count=count)

@app.route("/view_graphs")
def view():

    f5_path = request.args['f5_path']
    type = request.args['type']
    num = int(request.args['graph_num'])

    with open(f5_path+"/data_"+type+".tsv", 'rt') as data:
        for i, l in enumerate(data):
            if i == num:
                l = l.strip('\n')
                l = l.split('\t')
                fast5 = l[0]
                readID = l[1]
                if "." in l[4]:
                    sig = np.array([float(i) for i in l[4:]], dtype=float)
                else:
                    sig = np.array([int(i) for i in l[4:]], dtype=int)
    plt.ioff()
    graph = dict()
    html_graph = view_sig(sig, type, readID, fast5)
    graph['id'] = str(readID)
    graph['file'] = fast5
    graph['num'] = num
    #graph['json'] = json.dumps(dic)
    graph['html'] = Markup(html_graph)
    return render_template("view_graphs.html", f5_path=f5_path, type=type, graph=graph, count=i)

#@app.route("/test")
#def bkapp_page():
#    with pull_session(url="http://localhost:5006/sliders") as session:
#        session.document.root[0].children[1].title.text ="special sliders"
#        script = sever_session(session_id=session_id, url='special sliders')
#        return render_template("embed.html", script=script, template="Flask")

@app.route("/delete")
def delete():
    f5_path = request.args['f5_path']
    type = request.args['type']
    if os.path.isfile(f5_path+"/data_"+type+".tsv"):
        os.remove(f5_path+"/data_"+type+".tsv")
    return redirect("/")

def extract_f5_all(filename, type):
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
    multi = True
    raw = False
    if type == "raw":
        raw = True
    with h5py.File(filename, 'r') as hdf:
        
        reads = list(hdf.keys())
        if 'read' not in reads[1]:
            multi = False

        # single fast5 files
        if not multi:
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
                if not raw:
                    f5_dic['raw'] = np.array(f5_dic['raw'], dtype=int)
                    f5_dic['raw'] = convert_to_pA_numpy(f5_dic['raw'], digitisation, range, offset)
                    f5_dic['raw'] = np.round(f5_dic['raw'], 2)

            except:
                traceback.print_exc()
                sys.stderr.write("extract_fast5_all():failed to extract raw signal or fastq from {}".format(filename))
                f5_dic = {}

        # multi fast5 files
        else:
            for read in reads:
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
                    if not raw:
                        f5_dic[read]['raw'] = np.array(f5_dic[read]['raw'], dtype=int)
                        f5_dic[read]['raw'] = convert_to_pA_numpy(f5_dic[read]['raw'], digitisation, range, offset)
                        f5_dic[read]['raw'] = np.round(f5_dic[read]['raw'], 2)
                    
                except:
                    traceback.print_exc()
                    sys.stderr.write("extract_fast5_all():failed to read readID: {}".format(read))
    
    return f5_dic, multi

def convert_to_pA_numpy(d, digitisation, range, offset):
    raw_unit = range / digitisation
    return (d + offset) * raw_unit

#def scale_outliers(sig, args):
    ''' Scale outliers to within m stdevs of median '''
    ''' Remove outliers that don't fit within the specified bounds '''
#    k = (sig > args.lim_low) & (sig < args.lim_hi)
#    return sig[k]

def view_sig(sig, type, name, file):
    '''
    View the squiggle
    '''
    fig, ax = plt.subplots(figsize=(14,7))
    x = range(len(sig))
    y = sig.tolist()
    lines = ax.plot(x,y, marker="o", markersize=2)    

    plugins.connect(fig, plugins.PointLabelTooltip(lines[0],labels=y))
    #plt.autoscale()
    ax.set_xlabel("")
    
    if type == 'raw':
        ax.set_title("Raw signal for: {} | File: {}".format(name, file))
        ax.set_ylabel("Current - Not scaled")
    else:
        ax.set_title("Signal for: {} | File: {}".format(name, file))
        ax.set_ylabel("Current (pA)")       

    #graph = mpld3.fig_to_dict(fig)
    html_graph = mpld3.fig_to_html(fig)
    return html_graph
    

if __name__ == "__main__":
    app.run(port="8080", debug=True)
