#!/usr/bin/env python3

import sys
from Bio import SeqIO
from Bio import Align
from collections import OrderedDict
import csv
import os
import hashlib
from datetime import datetime



def get_ref2alt(ref_file, alt_file, out_prefix, save_psa):
    ### pairwise alignment: coord only (no mutation)
    ## currently doesn't handle multifasta entries and needs all sequences to be on the same strand!
    ## TODO: add validation os psa with minimum identity between sequences or quit!

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

    # TODO: '0' if before start of reference sequence and smpl_lentgh if after but with WARNING in both case !
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


def transform_gff(input_gff, coords_refGnm2smplGnm, out_prefix, include_metada_in_gff):
    with open(input_gff, 'r') as f_in_gff:
        lines = f_in_gff.readlines()

    new_gff_lines = []
    version_line = None

    # Catch gff-version line (if defined on first line of input file)
    if lines[0].startswith("##gff-version"):
        version_line = lines[0]
        lines = lines[1:]

    for line in lines:
        if line.startswith("##sequence-region"):
            parts = line.rstrip().split(" ")
            if len(parts) == 4:
                ref_chr, ref_start, ref_end = parts[1], int(parts[2]), int(parts[3])
                smpl_chr_start, smpl_pos_start = transform_individual_coordinate(ref_chr, ref_start, coords_refGnm2smplGnm)
                smpl_chr_end, smpl_pos_end = transform_individual_coordinate(ref_chr, ref_end, coords_refGnm2smplGnm)

                if (smpl_chr_start and smpl_chr_end and smpl_pos_start and smpl_pos_end) and (smpl_chr_start == smpl_chr_end):
                    new_sequence_region = f"##sequence-region {smpl_chr_start} {smpl_pos_start} {smpl_pos_end}\\n"
                    new_gff_lines.append(new_sequence_region)
                else:
                    print(f"Warning: Could not transform coordinates for sequence-region: {line.strip()}")
            continue

        if line.startswith("#"):
            # keeped in option?
            continue

        parts = line.rstrip().split('\t')

        if len(parts) > 1:
            if len(parts) < 9:
                print(f"Error: Non commented line with less than 9 columns: |{line.strip()}|")
                sys.exit(1)
            ref_chr, feature_start, feature_end = parts[0], int(parts[3]), int(parts[4])

            smpl_chr_start, smpl_pos_start = transform_individual_coordinate(ref_chr, feature_start, coords_refGnm2smplGnm)
            smpl_chr_end, smpl_pos_end = transform_individual_coordinate(ref_chr, feature_end, coords_refGnm2smplGnm)

            if not (smpl_chr_start and smpl_chr_end and smpl_pos_start and smpl_pos_end):
                print(f"Error: No correspondence found for feature location ('{line.strip()}')")
                sys.exit(1)
            elif smpl_chr_start != smpl_chr_end:
                print(f"Error: Feature spans multiple chromosomes ('{line.strip()}')")
                sys.exit(1)

            # Update feature with new positions
            parts[0] = smpl_chr_start
            parts[3] = str(smpl_pos_start)
            parts[4] = str(smpl_pos_end)

            new_gff_lines.append('\t'.join(parts) + '\\n')

    # Write the output GFF with metadata
    out_gff_file = out_prefix + "_transferredAnnotation.gff"
    with open(out_gff_file, "w") as f_out_gff:
        # Write the version line if present
        if version_line:
            f_out_gff.write(version_line)
        
        # Write metadata if requested
        if include_metada_in_gff:
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
            f_out_gff.write(f"# Metadata added by the script '{os.path.abspath(__file__)}'\\n")
            f_out_gff.write(f"#   -from: '{os.path.abspath(input_gff)}'\\n#          [hashlib.md5_4096:'{in_gff_hash}']\\n")
            f_out_gff.write(f"#   -by: '{os.path.abspath(__file__)}'\\n#        [hashlib.md5_4096:'{script_hash}']\\n")
            f_out_gff.write(f"#   -on: '{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}'\\n")

        # Write transformed GFF content
        f_out_gff.writelines(new_gff_lines)


def main():
    ref_fa="${annot_fa}"
    ref_gff="${annot_gff}"
    smpl_fa="${sample_fa}"
    out_prefix="${out_prefix}"
    #save_psa=\${save_psa}
    #include_metada_in_gff=\${include_metada_in_gff}
    save_psa=True
    include_metada_in_gff=True

    coords_refGnm2smplGnm = get_ref2alt(ref_fa, smpl_fa, out_prefix, save_psa)
    transform_gff(ref_gff, coords_refGnm2smplGnm, out_prefix, include_metada_in_gff)


if __name__ == "__main__":
    main()
