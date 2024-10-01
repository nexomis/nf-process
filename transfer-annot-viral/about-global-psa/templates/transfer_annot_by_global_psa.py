#!/usr/bin/env python3


from Bio import SeqIO
from Bio import Align
from collections import OrderedDict
import csv

from BCBio import GFF
from Bio.SeqFeature import FeatureLocation
import subprocess
import sys
import hashlib
from datetime import datetime
import os
import tempfile




def get_ref2alt(ref_file, alt_file, out_prefix, save_psa):
    ##### pairwise alignement: only coord (no mutation)
    #        - for now, not manage multifasta input and need thas all sequence are in same strand
    #          (add validation following minimal identity between sequences or stop !)

    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 5
    aligner.mismatch_score = -3
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    aligner.target_end_gap_score = -10
    aligner.query_end_gap_score = -10

    coords_refGnm2smplGnm = {}
    
    parsed_records = dict()
    parsed_records["ref"] = list(SeqIO.parse(ref_file, "fasta"))
    parsed_records["alt"] = list(SeqIO.parse(alt_file, "fasta"))
    records = dict()

    for key in parsed_records.keys():
        if len(parsed_records[key]) != 1:
            sys.exit("1 and only 1 contig per fasta file required")
        records[key] = parsed_records[key][0]

    alignments = aligner.align(records["ref"].seq, records["alt"].seq)
    best = alignments[0]
    ref_id = records["ref"].id
    alt_id = records["alt"].id
    print("aligning %s with %s" % (ref_id, alt_id))

    # Parse the alignment: coords
    alignment_array = best.__array__()
    ref = "".join(map(lambda x: x.decode("utf-8"), alignment_array[0]))
    alt = "".join(map(lambda x: x.decode("utf-8"), alignment_array[1]))
    
    if ref_id not in coords_refGnm2smplGnm:
        coords_refGnm2smplGnm[ref_id] = OrderedDict()
    # base1
    ref_pos = 1
    alt_pos = 1

    # todo: '0' if before start of reference sequence and smpl_lentgh if after but with WARNING in both case !
    # TODO: add -1 if out of ref (5' or 3') using alt/ref_pos and length of ref/alt + add del/ins after end to prevent gaps at end ?
    for i in range(best.length):
        ref_char = ref[i]
        alt_char = alt[i]
        if ref_char != '-' and alt_char != '-':
            coords_refGnm2smplGnm[ref_id][ref_pos] = (alt_id, alt_pos)
            ref_pos += 1
            alt_pos += 1
        elif ref_char == '-' and alt_char != '-':
            coords_refGnm2smplGnm[ref_id][ref_pos] = (alt_id, alt_pos)
            alt_pos += 1
        elif ref_char != '-' and alt_char == '-':
            coords_refGnm2smplGnm[ref_id][ref_pos] = (alt_id, alt_pos)
            ref_pos += 1

    # Save results on files
    refPos_to_smplPos_file = out_prefix + "_genomicCoords.csv"
    with open(refPos_to_smplPos_file, mode='w', newline='') as csv_file:
        fieldnames = ['refChr', 'refPos', 'smplChr', 'smplPos']
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        writer.writeheader()
        for ref_chr, positions in coords_refGnm2smplGnm.items():
            for ref_pos, (smpl_chr, smpl_pos) in positions.items():
                writer.writerow({
                    'refChr': ref_chr,
                    'refPos': ref_pos,
                    'smplChr': smpl_chr,
                    'smplPos': smpl_pos
                })
    
    if save_psa:
        alignment_file = out_prefix + "_" + ref_id + "_vs_" + alt_id + "_globalAlgn.txt"
        with open(alignment_file, 'w') as f_align:
            f_align.write("Best Alignment in plain text:\\n")
            f_align.write(best.format('clustal'))

    return coords_refGnm2smplGnm


def transform_individual_coordinate(ref_chr, ref_pos, coords_dict):
    if ref_chr in coords_dict and ref_pos in coords_dict[ref_chr]:
        return coords_dict[ref_chr][ref_pos]
    return None, None


### coords_refGnm2smplGnm: base1
### gff3: base1, start AND end inclusive
### BCBio.GFF: base0, start incluve but end exclusive 
### gff3.start = BCBio.GFF.start + 1 (base0 -> base1) | gff3.end = [BCBio.GFF.end + 1 (base0 -> base1) -1 (exlusive -> inclusive)]] = BCBio.GFF.end
def transform_gff(input_gff, coords_refGnm2smplGnm, out_prefix, include_metada_in_gff):
    with open(input_gff) as f_in_gff:
        gff_records = list(GFF.parse(f_in_gff))

    for record in gff_records:  # one record by chr
        smpl_chr = None
        for feature in record.features:
            #print(f"\\n>> record.id: '{record.id}' | feature.location.start+1: '{feature.location.start+1}' | feature.location.end: '{feature.location.end}'")
            smpl_chr_start, smpl_pos_start = transform_individual_coordinate(record.id, feature.location.start+1, coords_refGnm2smplGnm)
            smpl_chr_end, smpl_pos_end = transform_individual_coordinate(record.id, feature.location.end, coords_refGnm2smplGnm)
            #print(f"<<smpl_chr_start: '{smpl_chr_start}' | smpl_chr_end '{smpl_chr_end}' | smpl_pos_start '{smpl_pos_start}' | smpl_pos_end '{smpl_pos_end}'")
            if not (smpl_chr_start and smpl_chr_end and smpl_pos_start and smpl_pos_end):
                print(f"Error: no correspondance finding for feature location ('{feature}')")
                sys.exit(1)
            elif smpl_chr_start != smpl_chr_end:
                print(f"Error: multiple chromosome for feature (smpl_chr_start != smpl_chr_end: '{feature}')")
                sys.exit(1)
            else:
                if smpl_chr is None:
                    smpl_chr = smpl_chr_start
                elif smpl_chr != smpl_chr_start:
                    print(f"Error: multiple sample chromosome for one reference chromosome, not excepted for global alignment and recquire to check compatibility with gff transfer ('{record.id}' assigned to '{smpl_chr}' and to '{smpl_chr_start}')")
                    sys.exit(1)
                feature.location = FeatureLocation(smpl_pos_start, smpl_pos_end, feature.location.strand)
        record.id = smpl_chr
                
    # write output including metadata in case of 'include_metada_in_gff' True
    out_gff_file = out_prefix + "_transferedAnnotation.gff"
    if not include_metada_in_gff:
        with open(out_gff_file, "w") as f_out_gff:
            GFF.write(gff_records, f_out_gff)
    else:
        # tmpfile with raw out_gff without script metadata
        with tempfile.NamedTemporaryFile(delete=False, mode="w+") as temp_gff:
            GFF.write(gff_records, temp_gff)
            temp_gff_path = temp_gff.name
            
        with open(out_gff_file, "w") as f_out_gff:
            with open(temp_gff_path, "r") as temp_gff:
                # write metadata natively on temp_gff file
                while True:
                    line = temp_gff.readline()
                    if not line.startswith("#") or line.startswith("##sequence-region"):
                        break
                    f_out_gff.write(line)
                    
                # add metadata of script execution
                hash_md5 = hashlib.md5()
                with open(input_gff, "rb") as f:
                    for chunk in iter(lambda: f.read(4096), b""):
                        hash_md5.update(chunk)
                in_gff_hash = hash_md5.hexdigest()
                hash_md5 = hashlib.md5()
                with open(__file__, "rb") as f:
                    for chunk in iter(lambda: f.read(4096), b""):
                        hash_md5.update(chunk)
                script_hash = hash_md5.hexdigest()
                f_out_gff.write("# This GFF file has been reconstructed by transferring annotations (based on global pairwise alignment):\\n")
                f_out_gff.write(f"#    -from: '{os.path.abspath(input_gff)}'\\n#           [hashlib.md5_4096:'{in_gff_hash}']\\n")
                f_out_gff.write(f"#    -by:   '{__file__}'\\n#           [hashlib.md5_4096:'{script_hash}']\\n")
                f_out_gff.write(f"#    -on:   '{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}'\\n")
                
                # write the rest of temp_gff file
                f_out_gff.write(line)
                for line in temp_gff:
                    f_out_gff.write(line)
                    
        os.remove(temp_gff_path)



def main():
    ref_fa="${annot_fa}"
    ref_gff="${annot_gff}"
    smpl_fa="${sample_fa}"
    out_prefix="${meta.id}"
    #save_psa=\${save_psa}
    #include_metada_in_gff=\${include_metada_in_gff}
    save_psa=True
    include_metada_in_gff=True

    coords_refGnm2smplGnm = get_ref2alt(ref_fa, smpl_fa, out_prefix, save_psa)
    transform_gff(ref_gff, coords_refGnm2smplGnm, out_prefix, include_metada_in_gff)



if __name__ == "__main__":
    main()

