from __future__ import print_function
import argparse
import sys
import re
from statistics import mean
from collections import defaultdict

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile", help="A CosmicStructExport.tsv file", required=True, dest="infile")

    return parser.parse_args()

def make_dummy_header():
    h = '##fileformat=VCFv4.3\n'
    h += '##source=cosmic_structexport_to_vcf.py\n'
    h += '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV Type">\n'
    h += '##INFO=<ID=SVLEN,Number=1,Type=String,Description="SV length">\n'
    h += '##INFO=<ID=END,Number=1,Type=String,Description="END position of SV">\n'
    h += '##INFO=<ID=SPAN,Number=1,Type=String,Description="SPAN of SV">\n'
    h += '##INFO=<ID=CHR1_BP,Number=1,Type=String,Description="Breakpoint of SV on first CHROM">\n'
    h += '##INFO=<ID=CHR2_BP,Number=1,Type=String,Description="Breakpoint of SV on second CHROM">\n'
    h += '##INFO=<ID=CHR1,Number=1,Type=String,Description="First chromosome of SV">\n'
    h += '##INFO=<ID=CHR2,Number=1,Type=String,Description="Second chromosome of SV">\n'
    h += '##INFO=<ID=REF,Number=1,Type=String,Description="Reference genome of SV">\n'
    h += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"

    return h

def make_header_dict(line):
    h = defaultdict(int)
    r = defaultdict(str)
    index = 0
    for i in line.strip().split("\t"):
        h[i] = index
        r[index] = i
        index += 1

    return h, r

def info_to_string(info):
    s = []
    for i in info:
        if info[i] is not None and info[i] != "":
            s.append("=".join([i, str(info[i])]))

    return ";".join(s)

if __name__ == "__main__":

    args = parse_args()

    header_d = None
    reverse_d = None
    with open(args.infile, "r") as ifi:
        for line in ifi:
            if header_d is not None:
                tokens = line.strip().split("\t")
                mut_type = tokens[header_d["Mutation Type"]]
                
                svtype = None
                svtag = None
                svlen = None
                span = None
                chrom = None
                chrom2 = None
                pos = None
                pos2 = None
                end = None
                ref = None
                alt = None
                id_field = None
                qual = "99"
                filter_field = ""
                info = defaultdict(str)

                location_info = tokens[header_d["description"]]

                chrom = location_info.split(":")[0]
                location_info = ":".join(location_info.split(":")[1:])
                if mut_type  == "intrachromosomal with non-inverted orientation" or \
                        mut_type == "intrachromosomal with inverted orientation":
                    location_info = re.sub("_chr[0-9XYMT]*:", "_", location_info)
                if "_chr" in location_info:
                    chrom2 = re.findall("_chr[0-9XYMT]*:", location_info)[0].strip("_:")

                location_info = location_info.strip("g.o")
                stripped_location_info = location_info.strip("delinsinvdupbrkpttra")
                intervals = re.findall("\([0-9]*_[0-9XYMT]*\)", stripped_location_info)
                if len(intervals) == 0:
                    splits = stripped_location_info.split("_")
                    pos = splits[0]
                    end = splits[1]
                elif len(intervals) == 2:
                    split_intervals = [i.strip("()").split("_") for i in intervals]
                    mean_intervals = sorted([mean([int(j) for j in i]) for i in split_intervals])
                    pos = mean_intervals[0]
                    end = mean_intervals[1]
                else:
                    sys.stderr.write("Nonstandard start/length intervals: " + line + "\n")
                    continue
                   
                if "del" in location_info:
                    svtype = "DEL"
                    svtag = "<DEL>"
                   
                elif "ins" in location_info:
                    svtype = "INS"
                    svtag = "<INS>"
                    continue
                elif "inv" in location_info:
                    svtype = "INV"
                    svtag = "<INV>"

                elif "dup" in location_info:
                    svtype = "DUP"
                    svtag = "<DUP>"

                elif "bkpt" in location_info:
                        ## Requires special handling because for some reason this seemed
                        ## a reasonable format...
                        if "inverted" in mut_type:
                            svtype = "INV"
                            svtag = "<INV>"
                        elif not "inter" in mut_type.lower():
                            svtype = "DEL"
                            svtag = "<DEL>"
                        elif "inter" in mut_type.lower():
                            svtype = "TRA"
                            svtag = "<TRA>"

                elif "tra" in location_info:
                    continue
                else:
                    sys.stderr.write("Invalid SV type: " + location_info + "\n")
                    continue


                if not str(pos).isdigit() or not str(end).isdigit():
                    continue
 
                pos = int(pos)
                end = int(end)

                if pos > end:
                    tmp = pos
                    pos = end
                    end = tmp
                svlen = end - pos
                span = svlen
                if svtype == "DEL":
                    svlen = -1 * svlen

                info["SVLEN"] = svlen
                info["SVTYPE"] = svtype
                info["END"] = end
                info["SPAN"] = span
                if chrom2 is not None:
                    info["CHR2"] = chrom2
                info["MUTTYPE"] = str('"' + mut_type.replace(" ", "_") + '"')

                id_field = "."
                filter_field = "PASS"
                ref = "<N>"
                alt = svtag

                
                infos = info_to_string(info)
                if svtype is not None:
                    print("\t".join([str(i) for i in [chrom, pos, id_field, ref, alt, qual, filter_field, infos]]))

            else:
                header_d, reverse_d = make_header_dict(line)
                print(make_dummy_header())
