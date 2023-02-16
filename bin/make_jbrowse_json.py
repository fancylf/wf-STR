
import os
import json
import argparse


def check_file(path):
    if not os.path.exists(path):
        raise Exception(f"No Such File! {path}")
    else:
        return os.path.abspath(path)


def set_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--genome_fa')
    parser.add_argument('--location_chr')
    parser.add_argument('--location_start', type=int)
    parser.add_argument('--location_end', type=int)
    parser.add_argument('--bam')
    parser.add_argument('--template')

    return parser.parse_args()


def main():
    args = set_args()
    genome_fa = check_file(args.genome_fa)
    genome_fai = check_file(f"{genome_fa}.fai")
    bam = check_file(args.bam)
    bam_bai = check_file(f"{bam}.bai")
    location = "{}:{}..{}".format(args.location_chr, args.location_start, args.location_end)
    temp = check_file(args.template)

    config = json.loads(open(temp).read())
    config["assemblies"][0]["sequence"]["adapter"]["fastaLocation"]["uri"] = genome_fa
    config["assemblies"][0]["sequence"]["adapter"]["faiLocation"]["uri"] = genome_fai
    config["defaultSession"]["views"][0]["location"] = location
    config["tracks"][0]["adapter"]["bamLocation"]["uri"] = bam
    config["tracks"][0]["adapter"]["index"]["uri"] = bam_bai

    with open('jbrowse.json', 'w') as out:
        out.write(json.dumps(config, indent=4))


if __name__ == "__main__":
    main()
    
