
import argparse
import pandas as pd


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

    headers = ['#Name', 'Chromosome', 'Start position', 
               'End position', 'Motif', 'Reference repeats number', 
               'N', 'Region in genome', 'Gene', 'Pathogenic repeats number']
    df = pd.read_csv(info_file, sep='\t')
    db = pd.read_csv(database_file, header=None, names=headers)
    db.drop('N', axis=1, inplace=True)
    dfmerge = pd.merge(df, db, how='left', on='#Name')
    dfmerge = dfmerge.fillna('.')
    dfmerge['LP'] = dfmerge.apply(judge_positive, axis=1)
    dfmerge['STR score'] = ''

    return dfmerge


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--info', help="input GrandSTR info files")
    parser.add_argument('--pa', help="input pa csv file, 10 columns.")
    parser.add_argument('--outfile', help="output file. default: {info}.annot.xls")
    args = parser.parse_args()

    if args.outfile is None:
        args.outfile = f"{args.info}.annot.xls"

    df = add_grandstr_annot(args.info, args.pa)
    df.to_csv(args.outfile, sep='\t', index=False)


if __name__ == "__main__":
    main()
