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
        max_repeat = parse_repeat_string(ds['Details'])
        if max_repeat >= ds['Pathogenic repeats number']:
            return 'Yes'
        else:
            return 'No'

    df = pd.read_csv(info_file, sep='\t')
    db = pd.read_csv(database_file, sep='\t')
    dfmerge = pd.merge(df, db, how='left', on='#Name')
    dfmerge = dfmerge.fillna('.')
    dfmerge['LP'] = dfmerge.apply(judge_positive, axis=1)

    return dfmerge


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--info', help="GrandSTR info files")
    parser.add_argument('--db', default='grandstr_database.xls')
    args = parser.parse_args()

    df = add_grandstr_annot(args.info, args.db)
    df.to_csv(f"{args.info}.annot.xls", sep='\t', index=False)


if __name__ == "__main__":
    main()
