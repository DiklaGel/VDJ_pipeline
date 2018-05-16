##############################################################################
#                         Functions for use with                             #
# TraCeR - a tool to reconstruct TCR sequences from single-cell RNA-seq data #
#                                                                            #
# Please see README and LICENCE for details of use and licence conditions.   #
# This software was written by Mike Stubbington (ms31@sanger.ac.uk) from the #
# Teichmann Lab, EMBL-EBI and WTSI (www.teichlab.org). Latest versions are   #
# available for download at www.github.com/teichlab/tracer.                  #
#                                                                            #
#      Copyright (c) 2015, 2016 EMBL - European Bioinformatics Institute     #
#      Copyright (c) 2016 Genome Research Ltd.                               #
#      Author: M.J.T. Stubbington ms31@sanger.ac.uk                          #
##############################################################################

from __future__ import print_function

import copy
import os
import re
import subprocess
from collections import defaultdict

import Levenshtein
import six
from Bio import SeqIO
# import networkx as nx
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq

from gelSeqLib import core


def process_chunk(chunk):
    store_VDJ_rearrangement_summary = False
    store_junction_details = False
    store_alignment_summary = False
    store_hit_table = False
    store_sub_region = False
    alignment_summary = []
    hit_table = []
    looking_for_end = False
    return_dict = defaultdict(list)
    for line_x in chunk:
        if store_sub_region:
            sub_region_details = line_x.split("\t")
            for i in sub_region_details:
                return_dict['sub_region_details'].append(i)
            store_sub_region = False

        if store_VDJ_rearrangement_summary:
            VDJ_rearrangement_summary = line_x.split("\t")
            for i in VDJ_rearrangement_summary:
                return_dict['VDJ_rearrangement_summary'].append(i)
            store_VDJ_rearrangement_summary = False

        elif store_junction_details:
            junction_details = line_x.split("\t")
            for i in junction_details:
                return_dict["junction_details"].append(i)
            store_junction_details = False

        elif store_alignment_summary:
            if not line_x.startswith("#"):
                if line_x.startswith("Total"):
                    store_alignment_summary = False
                else:
                    return_dict['alignment_summary'].append(line_x)

        elif store_hit_table:
            if not looking_for_end:
                if not line_x.startswith("#"):
                    return_dict['hit_table'].append(line_x)
                    looking_for_end = True
            else:
                if line_x.startswith("#") or line_x.startswith("\n"):
                    store_hit_table = False
                else:
                    return_dict['hit_table'].append(line_x)

        elif line_x.startswith('# Query'):
            query_name = line_x.strip().split(" ")[2]
            #query_length = line_x.split(" ")[3]
            #return_dict['query_length'] = int(query_length.split("=")[1])
            return_dict['query_name'] = query_name

        elif line_x.startswith('# Sub-region sequence details'):
            store_sub_region = True

        elif line_x.startswith('# V-(D)-J rearrangement summary'):
            store_VDJ_rearrangement_summary = True

        elif line_x.startswith('# V-(D)-J junction details'):
            store_junction_details = True

        elif line_x.startswith('# Alignment summary'):
            store_alignment_summary = True

        elif line_x.startswith('# Hit table'):
            store_hit_table = True
    return (query_name, return_dict)


def find_possible_alignments(sample_dict, locus_names, cell_name, IMGT_seqs, output_dir, species, seq_method,
                             invariant_seqs, loci_for_segments, receptor, loci,
                             max_junc_string_length):
    alignment_dict = defaultdict(dict)
    recombinants = {}
    for locus in locus_names:
        recombinants[locus] = []

    # recombinants = {'TCRA': [], 'TCRB': []}
    for locus in locus_names:
        data_for_locus = sample_dict[locus]
        if data_for_locus is not None:
            for query_name, query_data in six.iteritems(data_for_locus):
                processed_hit_table = process_hit_table(query_name, query_data,
                                                        locus)

                if processed_hit_table is not None:
                    (returned_locus, good_hits,
                     rearrangement_summary) = processed_hit_table
                    junction_list = query_data['junction_details']

                    best_V = remove_allele_stars(
                        rearrangement_summary[0].split(",")[0])

                    junc_string = "".join(junction_list)
                    junc_string = remove_NA(junc_string)

                    locus_letter = returned_locus.split("_")[1]

                    if locus_letter in loci_for_segments['D']:
                        has_D = True
                    else:
                        has_D = False

                    if has_D:
                        best_J = remove_allele_stars(
                            rearrangement_summary[2].split(",")[0])
                    else:
                        best_J = remove_allele_stars(
                            rearrangement_summary[1].split(",")[0])

                    identifier = best_V + "_" + junc_string + "_" + best_J

                    # line attempting to add alignment summary to data for use
                    # with PCR comparisons
                    alignment_summary = query_data['alignment_summary']

                    all_V_names = [remove_allele_stars(x) for x in
                                   rearrangement_summary[0].split(',')]

                    if has_D:
                        all_J_names = [remove_allele_stars(x) for x in
                                       rearrangement_summary[2].split(',')]
                    else:
                        all_J_names = [remove_allele_stars(x) for x in
                                       rearrangement_summary[1].split(',')]

                    all_poss_identifiers = set()
                    for V in all_V_names:
                        for J in all_J_names:
                            i = V + "_" + junc_string + "_" + J
                            all_poss_identifiers.add(i)

                    # get original sequence from fasta file - needed for summary of reconstructed lengths.
                    # Only use the VDJ portion found by IgBLAST
                    fasta_file = "{output_dir}/{cell_name}_{locus}_collapsed.fasta".format(
                        output_dir=output_dir, cell_name=cell_name,locus=locus)
                    with open(fasta_file, 'rU') as fa:
                        for record in SeqIO.parse(fa, 'fasta'):
                            if query_name in record.id:
                                fasta_seq = record

                    if 'reversed' in good_hits[0][1]:
                        fasta_seq = fasta_seq.reverse_complement().seq
                    else:
                        fasta_seq = fasta_seq.seq
                    start_coord, end_coord = get_coords(good_hits)
                    fasta_seq = str(fasta_seq[start_coord:end_coord])

                    (imgt_reconstructed_seq, is_productive,
                     bestVJNames) = get_fasta_line_for_contig_imgt(
                        rearrangement_summary, junction_list, good_hits,
                        returned_locus, IMGT_seqs, cell_name,
                        query_name, species, loci_for_segments)
                    del (is_productive)
                    del (bestVJNames)

                    if seq_method == 'imgt':
                        (fasta_line_for_contig, is_productive,
                         bestVJNames) = get_fasta_line_for_contig_imgt(
                            rearrangement_summary, junction_list, good_hits,
                            returned_locus, IMGT_seqs, cell_name,
                            query_name, species, loci_for_segments)
                    cdr3 = data_for_locus[query_name]['sub_region_details']

                    if len(junc_string) < max_junc_string_length:
                        rec = core.Recombinant(contig_name=query_name,
                                          locus=returned_locus,
                                          identifier=identifier,
                                          all_poss_identifiers=all_poss_identifiers,
                                          productive=is_productive[0],
                                          stop_codon=is_productive[1],
                                          in_frame=is_productive[2], TPM=0.0,
                                          dna_seq=fasta_line_for_contig,
                                          hit_table=good_hits,
                                          summary=rearrangement_summary,
                                          junction_details=junction_list,
                                          best_VJ_names=bestVJNames,
                                          alignment_summary=alignment_summary,
                                          fasta_seq=fasta_seq,
                                          imgt_reconstructed_seq=imgt_reconstructed_seq,
                                          has_D=has_D,
                                               cdr3=cdr3)
                        recombinants[locus].append(rec)


    if recombinants:
        '''
        for locus, rs in six.iteritems(recombinants):
            # Adding code to collapse sequences with very low Levenshtein distances caused by confusion between
            # TRAVxD and TRAVx segments with different alignment lengths from
            # IgBlast.
            recombinants[locus] = collapse_close_sequences(rs, locus)
        '''
        cell = core.Cell(cell_name, recombinants, species=species, receptor=receptor,
                    loci=loci)
    else:
        cell = core.Cell(cell_name, None, species=species,
                    invariant_seqs=invariant_seqs, receptor=receptor, loci=loci)

    return (cell)


def get_fasta_line_for_contig_imgt(rearrangement_summary, junction_details,
                                   hit_table, locus, IMGT_seqs,
                                   sample_name, query_name, species,
                                   loci_for_segments):
    """
    Use first 258 bases of TRBC because they're the same between C1 and C2
    for TRGC use first 150 bases.
    Found by aligning the 4 C region transcripts and taking consensus.
    Ignored start of TCRG-C4-201 because it's only in that one.
    Use first 360 nt of TRBC1 because they're very nearly the same between
    TRBC1 and TRBCC2
    """

    found_best_V = False
    found_best_D = False
    found_best_J = False

    V_pattern = re.compile(r".+{potential_loci}V.+".format(
        potential_loci='[' + "".join(loci_for_segments['V']) + ']'))
    D_pattern = re.compile(r".+{potential_loci}D.+".format(
        potential_loci='[' + "".join(loci_for_segments['D']) + ']'))
    J_pattern = re.compile(r".+{potential_loci}J.+".format(
        potential_loci='[' + "".join(loci_for_segments['J']) + ']'))

    for hit in hit_table:
        segment = hit[2]
        V_match = V_pattern.search(segment)
        J_match = J_pattern.search(segment)
        if V_match and not found_best_V:
            # V_locus_key = "TR{}V".format(segment[2])
            V_locus_key = "_".join([locus, 'V'])
            best_V_name = segment
            # Remove forward slashes from shared A/D gene names to be the same as in the IMGT files.
            # segment = segment.replace("/", "_")
            best_V_seq = IMGT_seqs[V_locus_key][segment]

            # hit[11] is the end of the V sequence
            best_V_seq = best_V_seq[0:int(hit[11])]
            found_best_V = True
        elif J_match and not found_best_J:
            # J_locus_key = "TR{}J".format(segment[2])
            J_locus_key = "_".join([locus, 'J'])
            best_J_name = segment
            best_J_seq = IMGT_seqs[J_locus_key][segment]
            # hit 10 is the start of the J sequence
            best_J_seq = best_J_seq[int(hit[10]) - 1:]
            found_best_J = True

    junction = []

    parens_pattern = re.compile(r"\([CAGT]+\)")

    locus_letter = locus.split("_")[1]

    if locus_letter in loci_for_segments['D']:
        # junc_seqs = junction_details[1:3]
        VD_junc = junction_details[1]
        D_region = junction_details[2]
        DJ_junc = junction_details[3]
        if parens_pattern.search(VD_junc):
            VD_junc = re.sub(r'[\(\)]', '', VD_junc)
            length_in_parens = len(VD_junc)
            best_V_seq = best_V_seq[: -length_in_parens]
        if parens_pattern.search(DJ_junc):
            DJ_junc = re.sub(r'[\(\)]', '', DJ_junc)
            length_in_parens = len(DJ_junc)
            best_J_seq = best_J_seq[length_in_parens:]
        junc_seqs = [VD_junc, D_region, DJ_junc]

    else:
        VJ_junc = junction_details[1]
        # junctions in parentheses are represented in the coordinates of the matched segments.
        # Need to trim them then include the NTs in the junction
        if parens_pattern.search(VJ_junc):
            VJ_junc = re.sub(r'[\(\)]', '', VJ_junc)
            length_in_parens = len(VJ_junc)
            best_V_seq = best_V_seq[: -length_in_parens]
            best_J_seq = best_J_seq[length_in_parens:]
        junc_seqs = [VJ_junc]

    for seq in junc_seqs:
        seq = re.sub(r'[\(\)]', '', seq)
        if seq != "N/A":
            junction.append(seq)

    junction = "".join(junction)

    constant_seq = list(IMGT_seqs["_".join([locus, 'C'])].values())[0]

    # Editing IMGT V and J sequences to include any alterations from the
    # junction details
    V_end_seq = junction_details[0]
    J_start_seq = junction_details[-1]
    best_V_seq = best_V_seq[:-(len(V_end_seq))]

    best_V_seq = best_V_seq + V_end_seq
    best_J_seq = best_J_seq[len(J_start_seq):]
    best_J_seq = J_start_seq + best_J_seq

    full_rearrangement = best_V_seq + junction + best_J_seq + constant_seq
    productive_rearrangement = is_rearrangement_productive(
        best_V_seq + junction + best_J_seq + constant_seq[0:2])
    # fasta_line = ">chr={}__TCR{}_{}\n{}\n".format(sample_name, locus, query_name, full_rearrangement)

    bestVJ = [best_V_name, best_J_name]

    return (full_rearrangement, productive_rearrangement, bestVJ)


def is_rearrangement_productive(seq):
    # returns a tuple of three true/false values (productive, contains stop,
    # in-frame)
    seq_mod_3 = len(seq) % 3
    if seq_mod_3 == 0:
        in_frame = True
    else:
        in_frame = False

    seq = Seq(seq, IUPAC.unambiguous_dna)
    aa_seq = seq.translate()
    contains_stop = "*" in aa_seq

    if in_frame and not contains_stop:
        productive = True
    else:
        productive = False

    return (productive, contains_stop, in_frame)


def get_coords(hit_table):
    found_V = False
    found_J = False
    for entry in hit_table:
        if entry[0] == 'V':
            if not found_V:
                start = int(entry[8]) - 1
                found_V = True
        if entry[0] == 'J':
            if not found_J:
                end = int(entry[9])
                found_J = True
    return (start, end)


def remove_NA(junc_string):
    new_string = junc_string.replace("N/A", "")
    return (new_string)


def remove_allele_stars(segment):
    p = re.compile(r"(.+)\*\d+")
    m = p.search(segment)
    return (m.group(1))

def process_hit_table(query_name, query_data, locus):
    hit_table = query_data['hit_table']
    rearrangement_summary = query_data['VDJ_rearrangement_summary']

    e_value_cutoff = 5e-3

    found_V = set()
    found_D = set()
    found_J = set()

    good_hits = []

    segment_locus_pattern = re.compile(r"TR[AB]V.+DV.+")

    locus_name = locus.split("_")[1]

    for entry in hit_table:
        if not entry == "":
            entry = entry.split("\t")
            segment = entry[2]
            if segment_locus_pattern.search(segment):
                segment_locus = "AD"
            else:
                segment_locus = segment[2]
            segment_type = segment[3]
            e_value = float(entry[12])

            if locus_name in segment_locus:
                if e_value < e_value_cutoff:
                    if segment_type == "V":
                        found_V.add(locus)
                        good_hits.append(entry)
                    elif segment_type == "J":
                        found_J.add(locus)
                        good_hits.append(entry)
                else:
                    if segment_type == "D":
                        percent_identity = float(entry[3])
                        if percent_identity == 100:
                            found_D.add(locus)
                            good_hits.append(entry)

    if locus in found_V and locus in found_J:
        return (locus, good_hits, rearrangement_summary)
    else:
        return (None)



def is_rearrangement_productive(seq):
    # returns a tuple of three true/false values (productive, contains stop,
    # in-frame)
    seq_mod_3 = len(seq) % 3
    if seq_mod_3 == 0:
        in_frame = True
    else:
        in_frame = False

    seq = Seq(seq, IUPAC.unambiguous_dna)
    aa_seq = seq.translate()
    contains_stop = "*" in aa_seq

    if in_frame and not contains_stop:
        productive = True
    else:
        productive = False

    return (productive, contains_stop, in_frame)


def get_segment_name(name, pattern):
    match = pattern.search(name)
    number = match.group(1)
    if match.group(3):
        sub_number = match.group(3)
    else:
        sub_number = ""
    return (number)



def collapse_close_sequences(recombinants, locus):
    # pdb.set_trace()
    contig_names = [r.contig_name for r in recombinants]
    filtered_contig_names = [r.contig_name for r in recombinants]
    uncollapsible_contigs = []
    if len(recombinants) > 1:
        for i in range(len(recombinants) - 1):
            base_name = recombinants[i].contig_name
            base_seq = recombinants[i].imgt_reconstructed_seq
            base_V_segment = recombinants[i].best_VJ_names[0]
            base_J_segment = recombinants[i].best_VJ_names[1]

            base_id = recombinants[i].identifier
            base_junc = base_id.split("_")[1]
            base_e_value = float(recombinants[i].hit_table[0][-2])

            for j in range(i + 1, len(recombinants)):
                comp_name = recombinants[j].contig_name
                comp_seq = recombinants[j].imgt_reconstructed_seq
                comp_V_segment = recombinants[j].best_VJ_names[0]
                comp_J_segment = recombinants[j].best_VJ_names[1]
                comp_id = recombinants[j].identifier
                comp_junc = comp_id.split("_")[1]
                comp_e_value = float(recombinants[j].hit_table[0][-2])
                lev_dist = Levenshtein.distance(base_seq, comp_seq)
                # print("{}\t{}\t{}".format(base_id, comp_id, lev_dist))
                if lev_dist < 35 and not base_id == comp_id and base_name in filtered_contig_names \
                        and comp_name in filtered_contig_names:
                    # pdb.set_trace()
                    # define re pattern here to find TRAVx[DN] or TRDVx[DN]
                    # depending on locus
                    if locus == "TCR_A":
                        duplicate_pattern = re.compile(r"TRAV\d+[DN]")
                        segment_pattern = re.compile(
                            r"TRAV(\d+)([DN])?(-\d)?.+")
                        attempt_collapse = True
                    elif locus == "TCR_D":
                        duplicate_pattern = re.compile(r"DV\d+[DN]")
                        segment_pattern = re.compile(r"DV(\d+)([DN])?(-\d)?.+")
                        attempt_collapse = True
                    else:
                        uncollapsible_contigs.append(
                            "{}_vs_{}".format(base_name, comp_name))
                        attempt_collapse = False
                    if attempt_collapse and (
                        duplicate_pattern.search(
                            base_V_segment) or duplicate_pattern.search(
                            comp_V_segment)):
                        base_segment = get_segment_name(base_V_segment,
                                                        segment_pattern)
                        comp_segment = get_segment_name(comp_V_segment,
                                                        segment_pattern)
                        if base_segment == comp_segment:
                            # find alignment with lowest E value for V match
                            if base_e_value <= comp_e_value:
                                filtered_contig_names.remove(comp_name)
                            else:
                                filtered_contig_names.remove(base_name)
                        else:
                            uncollapsible_contigs.append(
                                "{}_vs_{}".format(base_name, comp_name))

                    else:
                        uncollapsible_contigs.append(
                            "{}_vs_{}".format(base_name, comp_name))

                elif lev_dist < 75 and not base_id == comp_id and base_name in filtered_contig_names \
                        and comp_name in filtered_contig_names:
                    if locus == "TCR_A":
                        duplicate_pattern = re.compile(r"TRAV\d+[DN]")
                        segment_pattern = re.compile(
                            r"TRAV(\d+)([DN])?(-\d)?.+")
                        attempt_collapse = True
                    elif locus == "TCR_D":
                        duplicate_pattern = re.compile(r"DV\d+[DN]")
                        segment_pattern = re.compile(r"DV(\d+)([DN])?(-\d)?.+")
                        attempt_collapse = True
                    else:
                        uncollapsible_contigs.append(
                            "{}_vs_{}".format(base_name, comp_name))
                        attempt_collapse = False
                    if attempt_collapse and (
                        duplicate_pattern.search(
                            base_V_segment) or duplicate_pattern.search(
                            comp_V_segment)):
                        base_segment = get_segment_name(base_V_segment,
                                                        segment_pattern)
                        comp_segment = get_segment_name(comp_V_segment,
                                                        segment_pattern)
                        if (base_segment == comp_segment) and (
                                base_junc == comp_junc) and (
                                base_J_segment == comp_J_segment):
                            # find alignment with lowest E value for V match
                            if base_e_value <= comp_e_value:
                                filtered_contig_names.remove(comp_name)
                            else:
                                filtered_contig_names.remove(base_name)
                        else:
                            uncollapsible_contigs.append(
                                "{}_vs_{}".format(base_name, comp_name))

                elif base_id == comp_id and base_name in filtered_contig_names and comp_name in filtered_contig_names:
                    if base_e_value <= comp_e_value:
                        filtered_contig_names.remove(comp_name)
                    else:
                        filtered_contig_names.remove(base_name)

    recombinants_to_delete = []

    for r in recombinants:
        if not r.contig_name in filtered_contig_names:
            recombinants_to_delete.append(r)
    print("recombinants_to_delete")
    print(str(recombinants_to_delete))
    [recombinants.remove(r) for r in recombinants_to_delete]

    return (recombinants)


'''
def make_cell_network_from_dna(cells, keep_unlinked, shape, dot, neato,
                               receptor, loci,
                               network_colours):
    G = nx.MultiGraph()
    # initialise all cells as nodes

    if shape == 'circle':
        for cell in cells:
            G.add_node(cell, shape=shape,
                       label=cell.html_style_label_for_circles(receptor, loci,
                                                               network_colours),
                       sep=0.4, fontname="helvetica neue")
    else:
        for cell in cells:
            G.add_node(cell, shape=shape,
                       label=cell.html_style_label_dna(receptor, loci,
                                                       network_colours),
                       fontname="helvetica neue")
    # make edges:
    for i in range(len(cells)):
        current_cell = cells[i]
        comparison_cells = cells[i + 1:]

        for locus in loci:
            col = network_colours[receptor][locus][0]

            # current_identifiers = current_cell.getMainRecombinantIdentifiersForLocus(locus)
            for comparison_cell in comparison_cells:
                shared_identifiers = 0
                if current_cell.recombinants[receptor][locus] is not None:
                    for current_recombinant in \
                            current_cell.recombinants[receptor][locus]:
                        current_id_set = current_recombinant.all_poss_identifiers
                        if comparison_cell.recombinants[receptor][
                                locus] is not None:
                            for comparison_recombinant in \
                                    comparison_cell.recombinants[receptor][locus]:
                                comparison_id_set = comparison_recombinant.all_poss_identifiers
                                if len(current_id_set.intersection(
                                        comparison_id_set)) > 0:
                                    shared_identifiers += 1

                if shared_identifiers > 0:
                    width = shared_identifiers * 2
                    G.add_edge(current_cell, comparison_cell, locus,
                               penwidth=width, color=col,
                               weight=shared_identifiers)

    deg = G.degree()

    to_remove = [n for n in deg if deg[n] == 0]

    if len(to_remove) < len(G.nodes()):
        if not shape == 'circle':
            G.remove_nodes_from(to_remove)
            drawing_tool = [dot, '-Gsplines=true', '-Goverlap=false',
                            '-Gsep=0.4']
        else:
            drawing_tool = [dot, '-Gsplines=true', '-Goverlap=false']
    else:
        drawing_tool = [neato, '-Gsplines=true', '-Goverlap=false']

    bgcolors = ['#8dd3c720', '#ffffb320', '#bebada20', '#fb807220', '#80b1d320',
                '#fdb46220', '#b3de6920', '#fccde520',
                '#d9d9d920', '#bc80bd20', '#ccebc520', '#ffed6f20']
    component_counter = 0
    component_groups = list()
    j = 0
    components = nx.connected_components(G)

    for component in components:
        members = list()
        if len(component) > 1:
            for cell in component:
                members.append(cell.name)
                G.node[cell]['style'] = 'filled'
                G.node[cell]['fillcolor'] = bgcolors[j]
                cell.bgcolor = bgcolors[j]

            if j < 11:
                j += 1
            else:
                component_counter += 1
                j = 0

        component_groups.append(members)

    return (G, drawing_tool, component_groups)


def draw_network_from_cells(cells, output_dir, output_format, dot, neato,
                            draw_graphs, receptor, loci, network_colours):
    cells = list(cells.values())
    network, draw_tool, component_groups = make_cell_network_from_dna(
        cells, False, "box", dot, neato, receptor, loci, network_colours)
    network_file = "{}/clonotype_network_with_identifiers.dot".format(
        output_dir)
    try:
        nx.write_dot(network, network_file)
    except AttributeError:
        import pydotplus
        nx.drawing.nx_pydot.write_dot(network, network_file)
    if draw_graphs:
        command = draw_tool + ['-o',
                               "{output_dir}/clonotype_network_with_identifiers.{output_format}".format(
                                   output_dir=output_dir,
                                   output_format=output_format), "-T",
                               output_format, network_file]
        subprocess.check_call(command)

    network, draw_tool, cgx = make_cell_network_from_dna(
        cells, False, "circle", dot, neato, receptor, loci, network_colours)
    network_file = "{}/clonotype_network_without_identifiers.dot".format(
        output_dir)
    try:
        nx.write_dot(network, network_file)
    except AttributeError:
        import pydotplus
        nx.drawing.nx_pydot.write_dot(network, network_file)
    if draw_graphs:
        command = draw_tool + ['-o',
                               "{output_dir}/clonotype_network_without_identifiers.{output_format}".format(
                                   output_dir=output_dir,
                                   output_format=output_format), "-T",
                               output_format, network_file]
        subprocess.check_call(command)
    return (component_groups)


def get_component_groups_sizes(cells, receptor, loci):
    cells = list(cells.values())
    G = nx.MultiGraph()
    # initialise all cells as nodes
    for cell in cells:
        G.add_node(cell)
    # make edges:
    for i in range(len(cells)):
        current_cell = cells[i]
        comparison_cells = cells[i + 1:]

        for locus in loci:

            for comparison_cell in comparison_cells:
                shared_identifiers = 0
                if current_cell.recombinants[receptor][locus] is not None:
                    for current_recombinant in \
                            current_cell.recombinants[receptor][locus]:
                        current_id_set = current_recombinant.all_poss_identifiers
                        if comparison_cell.recombinants[receptor][
                                locus] is not None:
                            for comparison_recombinant in \
                                    comparison_cell.recombinants[receptor][locus]:
                                comparison_id_set = comparison_recombinant.all_poss_identifiers
                                if len(current_id_set.intersection(
                                        comparison_id_set)) > 0:
                                    shared_identifiers += 1

                if shared_identifiers > 0:
                    width = shared_identifiers * 2
                    G.add_edge(current_cell, comparison_cell, locus,
                               penwidth=width, weight=shared_identifiers)

    deg = G.degree()

    to_remove = [n for n in deg if deg[n] == 0]

    # if len(to_remove) < len(G.nodes()):
    #    G.remove_nodes_from(to_remove)

    components = nx.connected_components(G)

    component_groups = list()

    singlets = []
    for component in components:
        members = list()
        if len(component) > 1:
            for cell in component:
                members.append(cell.name)
            component_groups.append(members)
        else:
            for cell in component:
                singlets.append(cell.name)

    clonotype_size_counter = Counter([len(x) for x in component_groups])
    clonotype_size_counter.update({1: len(singlets)})

    clonotype_sizes = []
    max_size = max(list(clonotype_size_counter.keys()))
    if max_size < 5:
        for x in range(1, max_size + 1):
            clonotype_sizes.append(clonotype_size_counter[x])
        zero_padding = 5 - len(clonotype_sizes)
        clonotype_sizes = clonotype_sizes + [0] * zero_padding
    else:
        for x in range(1, max_size + 1):
            clonotype_sizes.append(clonotype_size_counter[x])

    return (clonotype_sizes)

'''

def check_config_file(filename):
    if not os.path.isfile(filename):
        print()
        print("Couldn't find config file: {}".format(filename))
        print()
        exit(1)


def run_IgBlast(igblast, fasta, receptor, loci, output_dir, cell_name, species,index_location,
                ig_seqtype,aux_file_location,should_resume):
    print("##Running IgBLAST##")

    species_mapper = {
        'Mmus': 'mouse',
        'Hsap': 'human'
    }

    igblast_species = species_mapper[species]

    initial_locus_names = ["_".join([receptor, x]) for x in loci]
    locus_names = copy.copy(initial_locus_names)
    if should_resume:
        for locus in initial_locus_names:
            '''
            igblast_out = "{output_dir}/IgBLAST_output/{cell_name}_{receptor}_{locus}.IgBLASTOut".format(
                output_dir=output_dir, cell_name=cell_name,
                receptor=receptor, locus=locus)
            '''
            igblast_out = "{output_dir}/{cell_name}_{receptor}_{locus}.IgBLASTOut".format(
                output_dir=output_dir,
                receptor=receptor, locus=locus)
            if (os.path.isfile(igblast_out) and os.path.getsize(
                    igblast_out) > 0):
                locus_names.remove(locus)
                print(
                    "Resuming with existing IgBLAST output for {locus}".format(
                        locus=locus))

        if len(locus_names) == 0:
            return

    print("Performing IgBlast on {locus_names}".format(
        locus_names=locus_names))

    databases = {}
    for segment in ['V', 'D', 'J']:
        databases[segment] = "{}/{}_{}.fa".format(index_location, receptor,
                                                  segment)

    # Lines below suppress Igblast warning about not having an auxliary file.
    # Taken from
    # http://stackoverflow.com/questions/11269575/how-to-hide-output-of-subprocess-in-python-2-7
    DEVNULL = open(os.devnull, 'wb')

    for locus in locus_names:
        command = [igblast, '-germline_db_V', databases['V'],
                   '-germline_db_D', databases['D'],
                   '-germline_db_J', databases['J'], '-domain_system',
                   'imgt', '-organism', igblast_species,
                   '-ig_seqtype', ig_seqtype,'-auxiliary_data',aux_file_location, '-show_translation',
                   '-num_alignments_V', '2',
                   '-num_alignments_D', '2', '-num_alignments_J', '2',
                   '-outfmt', '7', '-query', fasta]
        igblast_out = "{output_dir}/{cell_name}_{locus}.IgBLASTOut".format(
                output_dir=output_dir, cell_name=cell_name,locus=locus)
        with open(igblast_out, 'w') as out:
            # print(" ").join(pipes.quote(s) for s in command)
            subprocess.check_call(command, stdout=out, stderr=DEVNULL)

    DEVNULL.close()

