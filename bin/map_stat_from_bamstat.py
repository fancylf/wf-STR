import argparse


def parse_samtools_bamstat(bamstat):
    d = {}

    for line in open(bamstat):
        if line.startswith('SN'):
            col = line.strip().split('\t')
            if col[1] == "raw total sequences:":
                d['total_reads'] = int(col[2])
            elif col[1] == "reads mapped:":
                d['map_reads'] = int(col[2])
            elif col[1] == "total length:":
                d['total_base'] = int(col[2])
            elif col[1] == "bases mapped (cigar):":
                d['map_base'] = int(col[2])
            else:
                pass

    return d


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample', help="sample name")
    parser.add_argument('--bamstat', help="samtools bamstat file")

    args = parser.parse_args()

    d = parse_samtools_bamstat(args.bamstat)

    headers = [
        "Sample_ID",
        "Pass_Reads",
        "Mapped_Reads",
        "Mapped_Reads_Rate(%)",
        "Pass_Bases",
        "Mapped_Bases",
        "Mapped_Bases_Rate(%)",
        "Depth(X)"
    ]

    with open(f"{args.sample}.map_stat.xls", 'w') as out:
        out.write('\t'.join(headers) + "\n")
        out.write("{}\t{:,}\t{:,}\t{:.2%}\t{:,}\t{:,}\t{:.2%}\t{}\n".format(
            args.sample,
            d['total_reads'],
            d['map_reads'],
            d['map_reads'] / d['total_reads'],
            d['total_base'],
            d['map_base'],
            d['map_base'] / d['total_base'],
            int(d['total_base'] / 3101804739)
        ))


if __name__ == "__main__":
    main()
    
