import pysam
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("--bam", nargs="+", dest="bams", help="BAM files")


def check_disjoint_cbs():
    global parser
    opts = parser.parse_args()
    barcodes = {}
    tot_alignments = 0

    for bam in opts.bams:
        print("reading " + bam)
        barcodes[bam] = {}
        with pysam.AlignmentFile(bam, "rb", check_sq=False) as input_alignments:
            for alignment in input_alignments:
                tot_alignments += 1
                if alignment.has_tag("CB"):
                    barcodes[bam][alignment.get_tag("CB")] = True

    for bam in opts.bams:
        print("checking " + bam)
        files = set(opts.bams)
        otherbams = files.difference(set([bam]))
        for cb in barcodes[bam].keys():
            for obam in otherbams:
                if cb in barcodes[obam]:
                    print("not a partition")
                    return

    print("tot alignments : ", tot_alignments)
    print("is a partition")
    return


if __name__ == "__main__":
    check_disjoint_cbs()
