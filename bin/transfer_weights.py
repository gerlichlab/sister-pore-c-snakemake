import os
import shutil
import argparse
import h5py

# define functions


def transferWeights(source, target, binsize):
    # get weights from source
    with h5py.File(source, 'r') as f:
        weights = f["resolutions"][binsize]["bins"]["weight"]
        weightArray = weights[:]
    # write weights into target
    with h5py.File(target, 'r+') as f:
        try:
            targetWeights = f["resolutions"][binsize]["bins"]["weight"]
            targetWeights[...] = weightArray
        except KeyError:
            f["resolutions"][binsize]["bins"]["weight"] = weightArray

# parse arguments

parser = argparse.ArgumentParser(description='Transfer weights from cis_and_trans to cis and trans')
parser.add_argument('--input_cis_and_trans')
parser.add_argument('--input_cis')
parser.add_argument('--input_trans')
parser.add_argument('--resolutions')
parser.add_argument('--output_cis_and_trans')
parser.add_argument('--output_cis')
parser.add_argument('--output_trans')
args = parser.parse_args()

# copy files to output -> needed for snakemake to work

shutil.copyfile(args.input_cis_and_trans, args.output_cis_and_trans)
shutil.copyfile(args.input_cis, args.output_cis)
shutil.copyfile(args.input_trans, args.output_trans)

parsed_resolutions = args.resolutions.split(",")

print(f"Transferring for resolutions: {parsed_resolutions}")

for binSize in parsed_resolutions:
    print(f"    Transferring from {args.input_cis_and_trans} to { args.output_cis} for binsize {binSize}")
    transferWeights(args.input_cis_and_trans, args.output_cis, binSize)
    print(f"    Transferring from {args.input_cis_and_trans} to { args.output_trans} for binsize {binSize}")
    transferWeights(args.input_cis_and_trans, args.output_trans, binSize)