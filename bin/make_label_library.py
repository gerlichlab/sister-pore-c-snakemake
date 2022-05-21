import argparse
import pickle
import pandas as pd

class BrdURead():
    def __init__(self, id, chrom, ref_start, ref_end, strand):
        self._read_id = id
        self._chrom = chrom
        self._ref_start = ref_start
        self._ref_end = ref_end
        self._strand = strand
        self._thymidines = []
    def add_thymidine(self, position_on_ref, probability):
        self._thymidines.append({"chrom": self._chrom,"pos": position_on_ref, "prob_brdu": probability})
    def get_id(self):
        return f"{self._read_id}_{self._chrom}_{self._ref_start}_{self._ref_end})"
    def get_length(self):
        return abs(self._ref_end - self._ref_start)
    def get_thymidines(self):
        return pd.DataFrame(self._thymidines)
    def is_labelled(self, prob_cutoff, brdu_cutoff):
        """Determines whether the read is labelled"""
        number_thy = len(self._thymidines)
        if number_thy == 0:
            return False
        number_labelled_thy = len([i for i in self._thymidines if i["prob_brdu"] > prob_cutoff])
        return (number_labelled_thy/number_thy) > brdu_cutoff
    def __repr__(self):
        return f"<BrdURead id:{self._read_id}>"


def build_alignment_index(detect_path, prob_cutoff=0.8, brdu_cutoff=0.1, limit=None):
    """reads in detect file and builds
    an index that holds information as to
    whether a read is labelled or not. Read_id + chrom + alignment position is used as key."""
    f = open(detect_path,'r')
    index = 0
    current_read = None
    alignment_index = dict()
    for line_number, line in enumerate(f):
        #ignore the header lines
        if line[0] == '#':
                continue
        #split the line into a list by whitespace
        splitLine = line.rstrip().split()
        if line[0] == '>':
            index += 1
            if current_read is not None:
                alignment_index[current_read.get_id()] = current_read.is_labelled(prob_cutoff, brdu_cutoff)
            if limit is not None and index > limit:
                return alignment_index
            readID = splitLine[0][1:]
            chromosome = splitLine[1]
            refStart = int(splitLine[2])
            refEnd = int(splitLine[3])
            strand = splitLine[4]
            current_read = BrdURead(
                    readID,
                    chromosome,
                    refStart,
                    refEnd,
                    strand
            )
        else:
            posOnRef = int(splitLine[0])
            probBrdU = float(splitLine[1])
            sixMerOnRef = splitLine[2]
            current_read.add_thymidine(posOnRef, probBrdU)
            #add these values to a container or do some processing here
    f.close()
    return alignment_index


# parse arguments

parser = argparse.ArgumentParser(description='Create label library containing mapping between read_id and labeling state')
parser.add_argument('--input')
parser.add_argument('--prob_cutoff', default=0.8)
parser.add_argument('--brdu_cutoff', default=0.1)
parser.add_argument('--output')
args = parser.parse_args()

# iterate through input and generate index

output_index = dict()

for input_detect in args.input:
    temp_index = build_alignment_index(input_detect, args.prob_cutoff, args.brdu_cutoff)
    output_index.update(temp_index)

# write output

with open(args.output, 'wb') as handle:
    pickle.dump(output_index, handle)