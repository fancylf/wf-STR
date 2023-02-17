#!/export/pipeline/RNASeq/Software/nbtools/bin/python

import os
import json
import argparse
from csv import DictReader
import plotly.express as px
import plotly
import pandas as pd
import jinja2


def mkdir(path):
    if not os.path.exists(path):
        os.makedirs(path)
    return os.path.abspath(path)

def parse_repeat_string(s):
    d = {}
    for i in s.split(','):
        k,v = i.split(':')
        d[int(k)] = int(v)
    x = list(range(1, max(d.keys())+11))
    y = list(map(lambda x: d[x] if x in d else 0, x))
    return x, y


def bar_plot(x, y, title, vline):
    df = pd.DataFrame({'repeats':x, 'reads number':y})
    fig = px.bar(df, x='repeats', y='reads number')
    fig.add_vline(x=vline, 
                  line_color='red', 
                  #annotation_text=str(vline), 
                  #annotation_position='bottom right',
                  line_dash='dot'
    )

    fig.update_layout(
        title_text=title,
        #title_x = 0.5,
        width=1000,
        height=600
    )

    return fig


def read_soft_version(txt):
    r = []
    for line in open(txt):
        soft, version = line.strip().split(',')
        r.append([soft, version])
    return r


def read_param_json(js):
    r = []
    for key,value in json.loads(open(js).read()).items():
        r.append([key, value])
    return r


def read_table_to_array(xls, header=True):
    fh = open(xls)
    if header:
        fh.readline()
    return [line.strip().split('\t') for line in fh]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--info', help="GrandSTR annotated info files")
    parser.add_argument('--mapstat', help="map stat xls")
    parser.add_argument('--outdir', help="output directory")
    parser.add_argument('--software', help="software versions txt")
    parser.add_argument('--param', help="parameter json file")
    parser.add_argument('--template', help="html template")
    args = parser.parse_args()

    #sample = os.path.basename(args.info).split('.')[0]
    
    for rec in DictReader(open(args.info), delimiter='\t'):
        if rec['LP'] == 'Yes':
            break
        if rec['#Name'] == 'STR046343':
            bak_rec = rec
    if rec['LP'] != 'Yes':
        rec = bak_rec

    x, y = parse_repeat_string(rec['Details'])
    title = f"STR_ID:{rec['#Name']} Gene:{rec['Gene']}"
    vline = int(rec['Pathogenic repeats number'])
    fig = bar_plot(x, y, title, vline)
    div = plotly.offline.plot(fig, output_type='div', include_plotlyjs=False)

    table_mapping = read_table_to_array(args.mapstat)
    table_soft = read_soft_version(args.software)
    table_param = read_param_json(args.param)
    temp = jinja2.Template(open(args.template).read())
    html = temp.render(
        divs=div, 
        table_soft=table_soft, 
        table_param=table_param,
        table_mapping = table_mapping
    )

    with open(os.path.join(args.outdir, 'STR_report.html'), 'w') as out:
        out.write(html)


if __name__ == "__main__":
    main()

