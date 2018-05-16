import os
import subprocess
from collections import defaultdict
from scipy.signal import argrelextrema
import itertools
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

def collapse_uniq_sequnces(fastq2):
    # reading all 15-mers in the current plate (current fastq2 - fastq file of read2), sorting them by frequency
    column = subprocess.getoutput(
        """gunzip -c %s | awk 'NR%s==2' | cut -b 1-15 | sort | uniq -c | sort -n """ % (fastq2, "%4")).split("\n")
    columns = [(int(column[i].strip().split(" ")[0]), column[i].strip().split(" ")[1][0:7],
                column[i].strip().split(" ")[1][7:15]) for i in range(0, len(column))]
    # creating a data frame from the 3-dim list above
    plate_seqs = pd.DataFrame(columns, columns=["num", "cell_barcode", "umi_barcode" ])
    high_confidence_barcodes = plate_seqs.sort_values(by= "num", ascending=False)

    return high_confidence_barcodes


def create_fasta_per_cell_stat(fastq1, fastq2, barcodes, output_dir):
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
            if ((barcodes["cell_barcode"] == cell_barcode) & (barcodes["umi_barcode"] == umi_barcode)).any():
                cell_name = str(barcodes[((barcodes["cell_barcode"] == cell_barcode) & (barcodes["umi_barcode"] == umi_barcode))].iloc[
                    0]["cell_name"])
                file_name = cell_name + "_" +cell_barcode + "_" + umi_barcode
                with open(file_name, 'a') as fa:
                    fasta_line = r1[line]
                    query_line = ">" + r1[line - 1][1:-1] + " " + cell_barcode + umi_barcode + "\n"
                    fa.write(query_line)
                    fa.write(fasta_line)
    os.remove(fastq1_dest)
    os.remove(fastq2_dest)


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


def create_fasta_per_cell(barcode_to_well, fastq1, fastq2, filtered_in_mapped_barcodes, output_dir):
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
                cell_name = barcode_to_well[cell_barcode]
                cell_fasta_file = output_dir + "/" + cell_name + ".fasta"
                with open(cell_fasta_file, 'a') as fa:
                    fasta_line = r1[line]
                    query_line = ">" + r1[line - 1][1:-1] + " " + cell_barcode + umi_barcode + "\n"
                    fa.write(query_line)
                    fa.write(fasta_line)
    os.remove(fastq1_dest)
    os.remove(fastq2_dest)


def split_to_cells(plate_name,wells_cells_file,output_dir,fastq1,fastq2,f):
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
    m["kmer_representative"] = m["cell_barcode"] + m["umi_barcode"]
    m["original_umi_barcode"] = ""
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
    log.write("filter kmers with less than " + str(t) + " repetitions" + "\n")
    un_mapped_barcodes = copy.deepcopy(un_mapped_barcodes.loc[un_mapped_barcodes["num"] >= t,])
    filtered_in_mapped_barcodes = copy.deepcopy(mapped_barcodes.loc[mapped_barcodes["num"] >= t,])
    log.write("total_reads_after_filtering\t" + str(filtered_in_mapped_barcodes["num"].sum()) + "\n")
    log.write("mapped_wells_filtered\t" + str(filtered_in_mapped_barcodes["Well_ID"].nunique()) + "\n")
    log.write("mapped_kmers_filtered\t" + str(len(mapped_barcodes)) + "\n")
    log.write("unmapped_kmers_filtered\t" + str(len(un_mapped_barcodes)) + "\n")
    log.write("total_unmapped_reads_after_filtering\t" + str(un_mapped_barcodes["num"].sum()) + "\n")


    sorted_wells = filtered_in_mapped_barcodes.groupby("Well_ID")["num"].max().sort_values(ascending=False).index
    for well in sorted_wells:
        well_table = filtered_in_mapped_barcodes[filtered_in_mapped_barcodes["cell_barcode"] == well_to_barcode[well]]
        for index, row in well_table.iterrows():
            unmapped = un_mapped_barcodes[un_mapped_barcodes["cell_barcode"].isin(barcode_dict[well_to_barcode[well]])]
            for index2, row2 in unmapped.iterrows():
                if hamming_distance(row["umi_barcode"], row2["umi_barcode"], 2) <= 1:
                    filtered_in_mapped_barcodes.loc[index2, ["Well_ID","well_coordinates"]] = row[["Well_ID","well_coordinates"]]
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
                if filtered_in_mapped_barcodes.loc[index,"original_umi_barcode"] == "":
                    dom_umi.append(row["umi_barcode"])
                    filtered_in_mapped_barcodes.loc[index,"original_umi_barcode"] = row["umi_barcode"]

    filtered_in_mapped_barcodes = filtered_in_mapped_barcodes.sort_index()
    # dropping out reads with similar/identical umi barcode (mapped to different wells) and similar umi barcode with lower frequency
    for field in ["original_umi_barcode","umi_barcode"]:
        for name,group in filtered_in_mapped_barcodes.groupby(by=field):
            if group["Well_ID"].nunique()>1:
                dom_well = group.iloc[0]["Well_ID"]
                filtered_in_mapped_barcodes = filtered_in_mapped_barcodes.drop(group[group["Well_ID"] != dom_well].index)
        log.write("total filterd mapped kmers after " + field + " filtering\t" + str(len(filtered_in_mapped_barcodes)) + "\n")
        log.write("total filterd mapped wells after " + field + " filtering\t" + str(filtered_in_mapped_barcodes["Well_ID"].nunique()) + "\n")


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


    log.write("total filterd mapped kmers after cell barcode similarity filtering\t" + str(len(filtered_in_mapped_barcodes)) + "\n")
    log.write("total filterd mapped wells after cell barcode similarity filtering\t" + str(
        filtered_in_mapped_barcodes["Well_ID"].nunique()) + "\n")

    filtered_in_mapped_barcodes[["Well_ID","well_coordinates","num","cell_barcode","umi_barcode","original_umi_barcode"]].to_csv(os.path.join(output_dir,"high_conf_before_control_cut.csv"),index=False)

    thresh = filtered_in_mapped_barcodes[filtered_in_mapped_barcodes["well_coordinates"].isin(["O1", "O2", "P1", "P2"])]["num"].max()
    if thresh > 0:
        log.write("control threshold\t" + str(thresh) +"\n")
        filtered_in_mapped_barcodes = filtered_in_mapped_barcodes[filtered_in_mapped_barcodes["num"] > thresh]

    filtered_in_mapped_barcodes[["Well_ID","well_coordinates","num","cell_barcode","umi_barcode","original_umi_barcode"]].to_csv(os.path.join(output_dir, "high_conf.csv"),index=False)
    log.write("total filterd mapped kmers after cutting the control threshold filtering\t" + str(len(filtered_in_mapped_barcodes)) +"\n")
    log.write("total filterd mapped wells after cutting the control threshold cell barcode similarity filtering\t" + str(
        filtered_in_mapped_barcodes["Well_ID"].nunique()) +"\n")

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
    log.close()
    create_fasta_per_cell(barcode_to_well, fastq1, fastq2, filtered_in_mapped_barcodes, output_dir)



"""""
    filt_in_wells = filtered_in_mapped_barcodes["Well_ID"].unique()
    filt_out_wells = [x for x in mapped_barcodes["Well_ID"].unique() if x not in filtered_in_mapped_barcodes["Well_ID"].unique()]
    filtered_out_mapped_barcodes = mapped_barcodes[(mapped_barcodes["num"] < t) & (mapped_barcodes["num"] > 4)]
    filtered_in_mapped_barcodes = mapped_barcodes[mapped_barcodes["num"] >= t]
    mapped_wells_filtered = len(pd.unique(filtered_in_mapped_barcodes["Well_ID"]))
    log.write("mapped_wells_filtered_after_first_collapse\t" + str(mapped_wells_filtered) + "\n")
    wells_out = pd.unique(filtered_out_mapped_barcodes["Well_ID"])
    filtered_out_wells = [x for x in wells_out if x not in pd.unique(filtered_in_mapped_barcodes["Well_ID"])]
    log.write("filtered_out_wells:\t" + str(len(filtered_out_wells)))
    # updating the filtered_out_mapped_barcodes
    changed_barcodes = dict()
    for index,row in filtered_out_mapped_barcodes.iterrows():
        i = find_relative(filtered_in_mapped_barcodes, row["kmer_representative"],3)
        if i!=-1:
            if filtered_in_mapped_barcodes.loc[i]["Well_ID"] != row["Well_ID"]:
                changed_barcodes[index] = i
                row["Well_ID"] = filtered_in_mapped_barcodes.loc[i]["Well_ID"]
                row["well_coordinates"] = filtered_in_mapped_barcodes.loc[i]["well_coordinates"]
                row["kmer_representative"] = filtered_in_mapped_barcodes.loc[i]["kmer_representative"]
                filtered_out_mapped_barcodes.loc[index] = row

    for index,row in filtered_out_mapped_barcodes.iterrows():
        i = find_umi_relative(filtered_out_mapped_barcodes[index:],row["umi_barcode"],row["cell_barcode"],1)
        if i!=-1:
            filtered_out_mapped_barcodes.loc[i]["kmer_representative"] = row["kmer_representative"]


    for name,group in filtered_out_mapped_barcodes.groupby("Well_ID"):
        for index,row in group.iterrows():
            if row["kmer_representative"] == row["cell_barcode"] + row["umi_barcode"]:


    kmer_distribution = filtered_mapped_barcodes.groupby("Well_ID")["num"].apply(list)
    kmer_distribution_2 = filtered_mapped_barcodes.groupby("Well_ID")["kmer_representative"].apply(list)
    pd.DataFrame({"kmer": kmer_distribution_2, "num": kmer_distribution}).to_csv(os.path.join(output_dir,"kmers_dist.csv"))


def new_split_by_cells(plate_name,wells_cells_file,output_dir,fastq1,fastq2,f):
    #fastq1 = "/home/labs/amit/weiner/Work/HANJAY/MiSeq/170621_M01759_0051_000000000-B8BWN/fastq/143_S7_L001_R1_001.fastq.gz"
    #fastq2 = "/home/labs/amit/weiner/Work/HANJAY/MiSeq/170621_M01759_0051_000000000-B8BWN/fastq/143_S7_L001_R2_001.fastq.gz"
    columns = [(int(column[i].strip().split(" ")[0]), column[i].strip().split(" ")[1][0:7],
                column[i].strip().split(" ")[1][7:15]) for i in range(0, len(column))]
    #map_cell_to_barcode = pd.read_csv("/home/labs/amit/diklag/files/wells_cells.txt", delimiter='\t', usecols=['Well_ID', 'well_coordinates', 'Cell_barcode', 'Amp_batch_ID'])
    map_cell_to_barcode = pd.read_csv(wells_cells_file, delimiter='\t',usecols = ['Well_ID','well_coordinates', 'Cell_barcode', 'Amp_batch_ID'])
    #map_cell_to_barcode = map_cell_to_barcode[map_cell_to_barcode['Amp_batch_ID'] == "AB2912"]
    map_cell_to_barcode = map_cell_to_barcode[map_cell_to_barcode['Amp_batch_ID'] == plate_name]
    plate_seqs = pd.DataFrame(columns, columns=["num", "cell_barcode", "umi_barcode"])
    x = np.array([int(column[i].strip().split(" ")[0]) for i in range(0, len(column))])
    # find the threshold number
    t = log_thresh(x)
    m = pd.merge(plate_seqs, map_cell_to_barcode, left_on="cell_barcode", right_on="Cell_barcode", how='outer')
    m = m.sort_values(by="num", ascending=False)
    mapped_barcodes = m[m["Well_ID"].notnull()]
    un_mapped_barcodes = m[m["Well_ID"].isnull()]
    un_mapped_barcodes = un_mapped_barcodes[un_mapped_barcodes["num"] > t]
    filtered_mapped_barcodes = mapped_barcodes
    filtered_mapped_barcodes["kmer_representative"] = filtered_mapped_barcodes["cell_barcode"] + filtered_mapped_barcodes["umi_barcode"]
    for index,row in un_mapped_barcodes.iterrows():
        cell_barcode = row["cell_barcode"]
        umi_barcode = row["umi_barcode"]
        kmer = cell_barcode+umi_barcode
        i = find_relative(filtered_mapped_barcodes,kmer)
        if i!=-1:
            new_row = row
            new_row["Well_ID"] = filtered_mapped_barcodes.loc[i]["Well_ID"]
            new_row["well_coordinates"] = filtered_mapped_barcodes.loc[i]["well_coordinates"]
            new_row["kmer_representative"] = filtered_mapped_barcodes.loc[i]["kmer_representative"]
            filtered_mapped_barcodes = filtered_mapped_barcodes.append(new_row)
    filtered_mapped_barcodes["kmer_total_num"] = filtered_mapped_barcodes["num"].groupby(
            filtered_mapped_barcodes["kmer_representative"]).transform('sum')
    filtered_mapped_barcodes = filtered_mapped_barcodes.sort_values(by="kmer_total_num",ascending=False)
    x = np.array(filtered_mapped_barcodes["kmer_total_num"].sort_values())
    hist, not_important = np.histogram(np.log(x), density=True, bins=int(np.ceil(np.log(x)[-1])))
    i = argrelextrema(hist, np.less)
    low_rate = filtered_mapped_barcodes[filtered_mapped_barcodes["kmer_total_num"] <= np.exp(i)]
    high_rate = filtered_mapped_barcodes[filtered_mapped_barcodes["kmer_total_num"] > np.exp(i)]
    for index,row in low_rate.iterrows():
        kmer = row["kmer_representative"]
        i = find_relative(high_rate,kmer)
        if i!=-1:
            row["kmer_representative"] = high_rate.loc[i]["kmer_representative"]
        high_rate = high_rate.append(row)
    high_rate["kmer_total_num"] = high_rate["num"].groupby(high_rate["kmer_representative"]).transform('sum')

"""""



