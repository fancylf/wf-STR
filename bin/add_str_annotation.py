import os
import argparse
import pandas as pd


BASEDIR = os.path.dirname(__file__)

def parse_repeat_string(s):
    d = {}
    for i in s.split(','):
        k,v = i.split(':')
        d[int(k)] = int(v)
    return max(d.keys())


def add_grandstr_annot(info_file, database_file):

    def judge_positive(ds):
        if ds['Pathogenic repeats number'] == '-':
            return '-'
        max_repeat = parse_repeat_string(ds['Details'])
        if max_repeat >= int(ds['Pathogenic repeats number']):
            return 'Yes'
        else:
            return 'No'

    df = pd.read_csv(info_file, sep='\t')
    db = pd.read_csv(database_file, sep='\t')
    dfmerge = pd.merge(df, db, how='left', on='#Name')
    dfmerge = dfmerge.fillna('.')
    dfmerge['LP'] = dfmerge.apply(judge_positive, axis=1)
    dfmerge['STR score'] = ''

    return dfmerge


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--info', help="GrandSTR info files")
    parser.add_argument('--db', help="database xls file, default: grandstr_database.xls")
    parser.add_argument('--outfile', help="output file. default: {info}.annot.xls")
    args = parser.parse_args()

    if args.db is None:
        args.db = os.path.join(BASEDIR, 'grandstr_database.xls')
    if args.outfile is None:
        args.outfile = f"{args.info}.annot.xls"

    df = add_grandstr_annot(args.info, args.db)
    df.to_csv(args.outfile, sep='\t', index=False)


if __name__ == "__main__":
    main()
