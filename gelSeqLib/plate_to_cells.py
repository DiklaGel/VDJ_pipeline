import os
import subprocess
from collections import defaultdict
from scipy.signal import argrelextrema
import itertools
from itertools import groupby
import numpy as np
import pandas as pd
import copy
import re
from Bio import pairwise2
from matplotlib import pyplot as plt
from gelSeqLib import io_func


################################################################
#  --------------------------to do-----------------------------#
# need to think of a better threshold function
# this function determines what is the minimal number for
################################################################

def threshold(groups,f):
    print("threshold= " + str(f))
    hist, not_important = np.histogram(groups, density=True, bins=groups[len(groups) - 1])
    sum = 0
    i = -1
    while sum + hist[i+1] < f:
        sum += hist[i+1]
        i += 1
    return i


def log_thresh(x):
    hist, not_important = np.histogram(np.log(x), density=True, bins=int(np.ceil(np.log(x)[-1])))
    i = argrelextrema(hist, np.less)[0][-1]
    return np.exp(i)

def gunzip_fastq(file,dest):
    subprocess.getoutput("""gunzip -c %s > %s  """ % (file, dest))

# simple hamming distance
def hamming_distance(str1,str2,max_dist):
    cost = 0
    for i in range(0,len(str1)):
        if str1[i] != str2[i]:
            cost += 1
            if cost >= max_dist:
                break
    return cost

def find_relative(df,kmer,distance):
    for i,row in df.iterrows():
        curr_kmer = row["cell_barcode"] + row["umi_barcode"]
        if hamming_distance(kmer,curr_kmer,distance+1) <=distance:
            return i
    return -1

def find_umi_relative(df,umi_barcode,cell_barcode,distance):
    for i,row in df.iterrows():
        if row["cell_barcode"] == cell_barcode:
            if hamming_distance(umi_barcode,row["umi_barcode"],distance+1) <=distance:
                return i
    return -1


def get_barcode_dict(map_cell_to_barcode,unmapped):
    well_to_barcode = pd.Series(map_cell_to_barcode.Cell_barcode.values, index=map_cell_to_barcode.Well_ID).to_dict()
    barcode_dict = defaultdict(list)
    b_list = list(map_cell_to_barcode["Cell_barcode"])
    for i in range(0, len(b_list)):
        barcode = b_list[i]
        for j in range(i + 1, len(b_list)):
            if hamming_distance(barcode, b_list[j], 3) <= 2:
                barcode_dict[barcode].append(b_list[j])
                barcode_dict[b_list[j]].append(barcode)
        for artificial_barcode in unmapped:
            if hamming_distance(barcode, artificial_barcode, 3) <= 2:
                barcode_dict[barcode].append(artificial_barcode)
    return barcode_dict, well_to_barcode


def create_fasta_per_cell_new(fastq1, fastq2, filtered_in_mapped_barcodes, output_dir,total_reads):
    fastq1_dest = os.path.join(output_dir, os.path.basename(fastq1).split(".gz")[0])
    gunzip_fastq(fastq1, fastq1_dest)
    fastq2_dest = os.path.join(output_dir, os.path.basename(fastq2).split(".gz")[0])
    gunzip_fastq(fastq2, fastq2_dest)
    final_output = pd.DataFrame(columns=["Well_ID", "cell_name",'plate_total_reads', "reads_freq", "umi_distribution", "unique_var_region"])
    with open(fastq1_dest) as f1, open(fastq2_dest) as f2:
        r1 = f1.readlines()
        #r2 = f2.readlines()
        for name,group in filtered_in_mapped_barcodes.groupby(by="well_coordinates"):
            cell_name = name
            cell_fasta_file = output_dir + "/" + cell_name + ".fasta"
            reads = []
            umi_dict = dict()
            for index,row in group.iterrows():
                cell_barcode = row["cell_barcode"]
                umi_barcode = row["umi_barcode"]
                original_umi_barcode = row["original_umi_barcode"]
                lines = subprocess.getoutput(
                """grep -n %s %s | cut -d ":" -f1 """ % (
                cell_barcode + umi_barcode,fastq2_dest))
                lines = lines.split("\n")
                lines = [int(i) - 1 for i in lines]
                umi_dict.update({key:original_umi_barcode for key in lines})
                reads.extend(lines)
            # generate a data frame of reads and their hyper variable region
            df2 = pd.DataFrame([(r1[i], r1[i][80:130],umi_dict[i]) for i in reads if len(r1[i]) > 130])
            if len(df2) == 0:
                continue
            df2.columns = ["read","hvr","original_umi_barcode"]
            df_full_reads = pd.DataFrame(df2.groupby(by=["hvr","original_umi_barcode", "read"]).size().sort_values(ascending=False))
            df_full_reads.columns = ["read_freq"]
            df_full_reads.reset_index(inplace=True)
            df_hvr = pd.DataFrame(df2.groupby(by=["hvr"]).size().sort_values(ascending=False).head())
            df_hvr.columns = ["hvr_freq"]
            df_hvr.reset_index(inplace=True)
            # peeking only reads with abundant hyper variable region
            m = pd.merge(df_full_reads, df_hvr, on="hvr", how='right').sort_values(by="read_freq",ascending=False).reset_index(drop=True)
            with open(cell_fasta_file, 'a') as fa:
                for index,row in m.iterrows():
                    query_line = ">" + str(index) + ":" + row["original_umi_barcode"] + "-" + str(row["read_freq"]) + "\n"
                    fa.write(query_line)
                    fa.write(row["read"])
            well_id = group["Well_ID"].iloc[0]
            final_output = final_output.append(
                    [{"Well_ID": well_id, "cell_name": cell_name, 'plate_total_reads': total_reads,
                      "reads_freq": m["read_freq"].sum(),
                      "umi_distribution": " ".join([str(int(count)) for count in
                                                    m.groupby("original_umi_barcode")["read_freq"].sum().sort_values(
                                                        ascending=False).values]),"unique_var_region": " ".join([str(int(count)) for count in
                                                    m.groupby("hvr")["read_freq"].sum().sort_values(
                                                        ascending=False).values])}])
            final_output = final_output.sort_values(by="reads_freq", ascending=False)
            final_output.to_csv(os.path.join(output_dir, "final_output.csv"), index=False)

    os.remove(fastq1_dest)
    os.remove(fastq2_dest)



def create_fasta_per_cell(fastq1, fastq2, filtered_in_mapped_barcodes, output_dir):
    fastq1_dest = os.path.join(output_dir, os.path.basename(fastq1).split(".gz")[0])
    gunzip_fastq(fastq1, fastq1_dest)
    fastq2_dest = os.path.join(output_dir, os.path.basename(fastq2).split(".gz")[0])
    gunzip_fastq(fastq2, fastq2_dest)
    with open(fastq1_dest) as f1, open(fastq2_dest) as f2:
        r1 = f1.readlines()
        r2 = f2.readlines()
        for line in range(1, len(r1), 4):
            cell_barcode = r2[line][0:7]
            umi_barcode = r2[line][7:15]
            if ((filtered_in_mapped_barcodes["cell_barcode"] == cell_barcode) & (
                        filtered_in_mapped_barcodes["umi_barcode"] == umi_barcode)).any():
                row = filtered_in_mapped_barcodes[(filtered_in_mapped_barcodes["cell_barcode"] == cell_barcode) & (
                        filtered_in_mapped_barcodes["umi_barcode"] == umi_barcode)].iloc[0]
                cell_name=row["well_coordinates"]
                cell_fasta_file = output_dir + "/" + cell_name + ".fasta"
                with open(cell_fasta_file, 'a') as fa:
                    fasta_line = r1[line]
                    query_line = ">" + r1[line - 1][1:-1] + " " + cell_barcode + umi_barcode + "\n"
                    fa.write(query_line)
                    fa.write(fasta_line)
    os.remove(fastq1_dest)
    os.remove(fastq2_dest)


def create_fasta_per_umi(fastq1, fastq2, filtered_in_mapped_barcodes, output_dir):
    fastq1_dest = os.path.join(output_dir, os.path.basename(fastq1).split(".gz")[0])
    gunzip_fastq(fastq1, fastq1_dest)
    fastq2_dest = os.path.join(output_dir, os.path.basename(fastq2).split(".gz")[0])
    gunzip_fastq(fastq2, fastq2_dest)
    with open(fastq1_dest) as f1, open(fastq2_dest) as f2:
        r1 = f1.readlines()
        r2 = f2.readlines()
        for line in range(1, len(r1), 4):
            cell_barcode = r2[line][0:7]
            umi_barcode = r2[line][7:15]
            if ((filtered_in_mapped_barcodes["cell_barcode"] == cell_barcode) & (
                        filtered_in_mapped_barcodes["umi_barcode"] == umi_barcode)).any():
                row = filtered_in_mapped_barcodes[(filtered_in_mapped_barcodes["cell_barcode"] == cell_barcode) & (
                        filtered_in_mapped_barcodes["umi_barcode"] == umi_barcode)].iloc[0]
                original_umi_barcode = row["original_umi_barcode"]
                cell_name=row["well_coordinates"]
                fasta_name = cell_name + "-" + original_umi_barcode
                cell_fasta_file = output_dir + "/" + fasta_name + ".fasta"
                with open(cell_fasta_file, 'a') as fa:
                    fasta_line = r1[line]
                    query_line = ">" + r1[line - 1][1:-1] + " " + cell_barcode + umi_barcode + "\n"
                    fa.write(query_line)
                    fa.write(fasta_line)
    os.remove(fastq1_dest)
    os.remove(fastq2_dest)


def create_fasta_per_umi_and_cell(fastq1, fastq2, filtered_in_mapped_barcodes, output_dir):
    fastq1_dest = os.path.join(output_dir, os.path.basename(fastq1).split(".gz")[0])
    gunzip_fastq(fastq1, fastq1_dest)
    fastq2_dest = os.path.join(output_dir, os.path.basename(fastq2).split(".gz")[0])
    gunzip_fastq(fastq2, fastq2_dest)
    cells_dir = os.path.join(output_dir,"cells")
    umis_dir = os.path.join(output_dir,"cells_umis")
    io_func.makeOutputDir(cells_dir)
    io_func.makeOutputDir(umis_dir)
    with open(fastq1_dest) as f1, open(fastq2_dest) as f2:
        r1 = f1.readlines()
        r2 = f2.readlines()
        for line in range(1, len(r1), 4):
            cell_barcode = r2[line][0:7]
            umi_barcode = r2[line][7:15]
            if ((filtered_in_mapped_barcodes["cell_barcode"] == cell_barcode) & (
                        filtered_in_mapped_barcodes["umi_barcode"] == umi_barcode)).any():
                row = filtered_in_mapped_barcodes[(filtered_in_mapped_barcodes["cell_barcode"] == cell_barcode) & (
                        filtered_in_mapped_barcodes["umi_barcode"] == umi_barcode)].iloc[0]
                original_umi_barcode = row["original_umi_barcode"]
                cell_name=row["well_coordinates"]
                write_fasta_line(cell_barcode, cell_name, r1[line], r1[line-1],cells_dir, umi_barcode)
                name = cell_name + "-" + original_umi_barcode
                write_fasta_line(cell_barcode,name,r1[line],r1[line-1], umis_dir, umi_barcode)
    os.remove(fastq1_dest)
    os.remove(fastq2_dest)

def write_fasta_line(cell_barcode, fasta_name, line,prev_line ,output_dir, umi_barcode):
    cell_fasta_file = output_dir + "/" + fasta_name + ".fasta"
    with open(cell_fasta_file, 'a') as fa:
        fasta_line = line
        query_line = ">" + prev_line[1:-1] + " " + cell_barcode + umi_barcode + "\n"
        fa.write(query_line)
        fa.write(fasta_line)


def split_to_cells(plate_name,wells_cells_file,output_dir,fastq1,fastq2,f):
    ## In order to reduce the noise and running time, we need to drop out unwanted reads
    # first way to handle it: filter by read2 (by kmers):
    ## 1. drop out reads with kmers frequency less than 4
    ## 2. droup out reads with same (similar or identitcal) umi barcode but different cell barcode and keep only the reads with the abundant cell barcode
    # second way to handle it: filter by read1 (by the gene sequene)
    ##

    log_file = os.path.join(output_dir,"split_log.log")
    log = open(log_file, 'w')
    column = subprocess.getoutput(
        """gunzip -c %s | awk 'NR%s==2' | cut -b 1-15 | sort | uniq -c | sort -n """ % (fastq2, "%4")).split("\n")
    columns = [(int(column[i].strip().split(" ")[0]), column[i].strip().split(" ")[1][0:7],
                column[i].strip().split(" ")[1][7:15]) for i in range(0, len(column))]
    plate_seqs = pd.DataFrame(columns, columns=["num", "cell_barcode", "umi_barcode"])
    total_reads = plate_seqs["num"].sum()
    log.write("total_reads\t" + str(total_reads) + "\n")
    #t = plate_seqs["num"].quantile(f)
    t = 4
    # reading wells cells file (mapping from cell barcode to well id/well coordinates
    map_cell_to_barcode = pd.read_csv(wells_cells_file, delimiter='\t',usecols = ['Well_ID','well_coordinates', 'Cell_barcode', 'Amp_batch_ID'])
    map_cell_to_barcode = map_cell_to_barcode.loc[map_cell_to_barcode['Amp_batch_ID'] == plate_name,]

    m = pd.merge(plate_seqs, map_cell_to_barcode, left_on="cell_barcode", right_on="Cell_barcode", how='outer')
    m = m.sort_values(by="num", ascending=False)
    m = m.reset_index(drop=True)
    m["original_umi_barcode"] = None
    mapped_barcodes = copy.deepcopy(m.loc[m["Well_ID"].notnull(),])
    un_mapped_barcodes = copy.deepcopy(m.loc[m["Well_ID"].isnull(),])
    # generating barcode to well dict
    barcode_dict, well_to_barcode = get_barcode_dict(map_cell_to_barcode,un_mapped_barcodes["cell_barcode"].unique())# generate a merged table
    log.write("mapped_reads\t" + str(mapped_barcodes["num"].sum()) + "\n")
    log.write("un_mapped_reads\t" + str(un_mapped_barcodes["num"].sum()) + "\n")
    quantile = ",".join([str(x) + ":" + str(int(m["num"].quantile(x))) for x in np.linspace(0.9, 0.99, 10)])
    log.write("quantile\t" + quantile + "\n")
    log.write("mapped_wells_not_filtered\t" + str(mapped_barcodes["Well_ID"].nunique()) + "\n")
    log.write("mapped_kmers_not_filtered\t" + str(len(mapped_barcodes)) + "\n")
    log.write("unmapped_kmers_not_filtered\t" + str(len(un_mapped_barcodes)) + "\n")
    thresh = mapped_barcodes[mapped_barcodes["well_coordinates"].isin(["O1", "O2", "P1", "P2"])]["num"].max()
    if thresh > 0:
        log.write("control_threshold_not_filtered\t" + str(thresh) +"\n")
    else:
        log.write("control_threshold_not_filtered\t" + str(0) + "\n")

    log.write("filter_kmers_with_less_than " + str(t) + " repetitions" + "\n")
    un_mapped_barcodes = copy.deepcopy(un_mapped_barcodes.loc[un_mapped_barcodes["num"] >= t,])
    filtered_in_mapped_barcodes = copy.deepcopy(mapped_barcodes.loc[mapped_barcodes["num"] >= t,])
    log.write("total_reads_after_filtering\t" + str(filtered_in_mapped_barcodes["num"].sum()) + "\n")
    log.write("mapped_wells_filtered\t" + str(filtered_in_mapped_barcodes["Well_ID"].nunique()) + "\n")
    log.write("mapped_kmers_filtered\t" + str(len(mapped_barcodes)) + "\n")
    log.write("unmapped_kmers_filtered\t" + str(len(un_mapped_barcodes)) + "\n")
    log.write("total_unmapped_reads_after_filtering\t" + str(un_mapped_barcodes["num"].sum()) + "\n")

    # iterating on unmapped reads with cell barcode similar to mapped one
    sorted_wells = filtered_in_mapped_barcodes.groupby("Well_ID")["num"].max().sort_values(ascending=False).index
    for well in sorted_wells:
        well_table = filtered_in_mapped_barcodes[filtered_in_mapped_barcodes["cell_barcode"] == well_to_barcode[well]]
        for index, row in well_table.iterrows():
            unmapped = un_mapped_barcodes[un_mapped_barcodes["cell_barcode"].isin(barcode_dict[well_to_barcode[well]])]
            for index2, row2 in unmapped.iterrows():
                if hamming_distance(row["umi_barcode"], row2["umi_barcode"], 2) <= 1:
                    filtered_in_mapped_barcodes.loc[index2, ["Well_ID","well_coordinates","original_umi_barcode"]] = row[["Well_ID","well_coordinates","original_umi_barcode"]]
                    filtered_in_mapped_barcodes.loc[index2, ["num","cell_barcode","umi_barcode"]] = row2[["num","cell_barcode","umi_barcode"]]
                    un_mapped_barcodes = un_mapped_barcodes.drop(index2)
    filtered_in_mapped_barcodes[["Well_ID","well_coordinates","num","cell_barcode","umi_barcode","original_umi_barcode"]].to_csv(os.path.join(output_dir,"high_conf_before_filtering.csv"),index=False)
    log.write("total_reads_after_adding_unmapped_barcodes\t" + str(filtered_in_mapped_barcodes["num"].sum()) + "\n")
    log.write("mapped_wells_after_adding_unmapped_barcodes\t" + str(filtered_in_mapped_barcodes["Well_ID"].nunique()) + "\n")
    log.write("mapped_kmers_after_adding_unmapped_barcodes\t" + str(len(mapped_barcodes)) + "\n")
    log.write("unmapped_kmers_after_adding_unmapped_barcodes\t" + str(len(un_mapped_barcodes)) + "\n")
    log.write("total_unmapped_reads_after_adding_unmapped_barcodes\t" + str(un_mapped_barcodes["num"].sum()) + "\n")
    # for each umi in each well - decide what is the original sequence
    filtered_in_mapped_barcodes = filtered_in_mapped_barcodes.sort_index()
    for name, group in filtered_in_mapped_barcodes.groupby(by="Well_ID"):
        if group["umi_barcode"].nunique()>1:
            dom_umi = [group.iloc[0]["umi_barcode"]]
            for index,row in group.iterrows():
                for dom in dom_umi:
                    if hamming_distance(row["umi_barcode"],dom,3) <=2:
                        filtered_in_mapped_barcodes.loc[index,"original_umi_barcode"] = dom
                        break
                if filtered_in_mapped_barcodes.loc[index,"original_umi_barcode"] is None:
                    dom_umi.append(row["umi_barcode"])
                    filtered_in_mapped_barcodes.loc[index,"original_umi_barcode"] = row["umi_barcode"]
        else:
            for index, row in group.iterrows():
                filtered_in_mapped_barcodes.loc[index,"original_umi_barcode"] = row["umi_barcode"]


    filtered_in_mapped_barcodes = filtered_in_mapped_barcodes.sort_index()
    # dropping out reads with similar/identical umi barcode (mapped to different wells) and similar umi barcode with lower frequency
    for field in ["original_umi_barcode","umi_barcode"]:
        for name,group in filtered_in_mapped_barcodes.groupby(by=field):
            if group["Well_ID"].nunique()>1:
                dom_well = group.iloc[0]["Well_ID"]
                filtered_in_mapped_barcodes = filtered_in_mapped_barcodes.drop(group[group["Well_ID"] != dom_well].index)
        log.write("total_filterd_mapped_kmers_after_" + field + "_filtering\t" + str(len(filtered_in_mapped_barcodes)) + "\n")
        log.write("total_filterd_mapped_wells_after_" + field + "_filtering\t" + str(filtered_in_mapped_barcodes["Well_ID"].nunique()) + "\n")

    # dropping out reads with similar cell barcode (mapped to different wells) and similar umi barcode with lower frequency
    sorted_wells = filtered_in_mapped_barcodes.groupby("Well_ID")["num"].max().sort_values(ascending=False).index
    for well in sorted_wells:
        well_table = filtered_in_mapped_barcodes[filtered_in_mapped_barcodes["cell_barcode"] == well_to_barcode[well]]
        for index, row in well_table.iterrows():
            sub_table = filtered_in_mapped_barcodes[(
                filtered_in_mapped_barcodes["cell_barcode"].isin(barcode_dict[well_to_barcode[well]])) & (filtered_in_mapped_barcodes["Well_ID"] != well)]
            for index2, row2 in sub_table.iterrows():
                if hamming_distance(row["umi_barcode"],row2["umi_barcode"],2) <=1:
                    if row["num"] >= row2["num"]:
                        filtered_in_mapped_barcodes = filtered_in_mapped_barcodes.drop(index2)
                    else:
                        filtered_in_mapped_barcodes = filtered_in_mapped_barcodes.drop(index)
                        break

    log.write("total_filterd_mapped_kmers_after_cell_barcode_similarity_filtering\t" + str(len(filtered_in_mapped_barcodes)) + "\n")
    log.write("total_filterd_mapped_wells_after_cell_barcode_similarity_filtering\t" + str(
        filtered_in_mapped_barcodes["Well_ID"].nunique()) + "\n")

    filtered_in_mapped_barcodes[["Well_ID","well_coordinates","num","cell_barcode","umi_barcode","original_umi_barcode"]].to_csv(os.path.join(output_dir,"high_conf_before_control_cut.csv"),index=False)

    old_thresh = 0 if thresh <= 0 else thresh
    thresh = filtered_in_mapped_barcodes[filtered_in_mapped_barcodes["well_coordinates"].isin(["O1", "O2", "P1", "P2"])]["num"].max()
    if thresh > 0:
        log.write("control_threshold_after_filtering\t" + str(thresh) +"\n")
    else:
        log.write("control_threshold_after_filtering\t" + str(old_thresh) + "\n")

    """
    filtered_in_mapped_barcodes = filtered_in_mapped_barcodes[filtered_in_mapped_barcodes["num"] > thresh]
    filtered_in_mapped_barcodes[["Well_ID","well_coordinates","num","cell_barcode","umi_barcode","original_umi_barcode"]].to_csv(os.path.join(output_dir, "high_conf.csv"),index=False)
    log.write("total_filterd_mapped_kmers_after_cutting_the_control_threshold_filtering\t" + str(len(filtered_in_mapped_barcodes)) +"\n")
    log.write("total_filterd_mapped_wells_after_cutting_the_control_threshold_cell_barcode_similarity_filtering\t" + str(
        filtered_in_mapped_barcodes["Well_ID"].nunique()) +"\n")
    """
    """
    cells_dir = os.path.join(output_dir,"cells")
    umis_dir = os.path.join(output_dir,"cells_umis")
    io_func.makeOutputDir(cells_dir)
    io_func.makeOutputDir(umis_dir)
    """

    map_cell_to_barcode = filtered_in_mapped_barcodes[["cell_barcode", "Well_ID","well_coordinates"]].drop_duplicates()
    barcode_to_well = pd.Series(map_cell_to_barcode.well_coordinates.values,index=map_cell_to_barcode.cell_barcode).to_dict()

    final_output = pd.DataFrame(columns=["Well_ID","cell_name","reads_freq","umi_distribution"])
    groups = filtered_in_mapped_barcodes.groupby("well_coordinates")
    for cell_name, cell_group in groups:
        well_id = cell_group.iloc[0]["Well_ID"]
        final_output = final_output.append([{"Well_ID": well_id, "cell_name": cell_name, 'plate_total_reads': total_reads,
                      "reads_freq": cell_group["num"].sum(),
                      "umi_distribution":" ".join([str(int(count)) for count in cell_group.groupby("original_umi_barcode")["num"].sum().sort_values(ascending=False).values])}])
    final_output = final_output.sort_values(by="reads_freq",ascending=False)
    final_output.to_csv(os.path.join(output_dir,"final_output.csv"),index = False)

    """
    final_output = pd.DataFrame(columns=["Well_ID","cell_name","reads_freq","umi_distribution"])
    groups = filtered_in_mapped_barcodes.groupby("well_coordinates")
    for cell_name, cell_group in groups:
        well_id = cell_group.iloc[0]["Well_ID"]
        for umi_name, umi_group in cell_group.groupby("original_umi_barcode"):
            final_output = final_output.append([{"Well_ID": well_id,"cell_name":cell_name+"-"+umi_name, 'plate_total_reads': total_reads,
                          "reads_freq": umi_group["num"].sum(),
                                                 "umi_distribution":" ".join([row["cell_barcode"]+row["umi_barcode"]+":"+str(int(row["num"]))
                                                                              for i,row in umi_group.sort_values(by="num",ascending=False).iterrows()])}])
    final_output = final_output.sort_values(by="reads_freq",ascending=False)
    final_output.to_csv(os.path.join(umis_dir,"final_output.csv"),index = False)
    """

    log.close()

    create_fasta_per_cell_new(fastq1, fastq2, filtered_in_mapped_barcodes, output_dir,total_reads)

