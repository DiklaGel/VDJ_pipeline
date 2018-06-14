import os
import subprocess
from collections import defaultdict, Counter
import numpy as np
import pandas as pd
import copy


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

def consensus(reads,freqs):
    bestseqs = [[]]
    #c = []
    max_len = Counter([len(c) for c in list(reads)]).most_common(1)[0][0]
    for i in range(0,max_len):
        c = Counter({N:sum([freqs[j] for j in range(0,len(freqs)) if len(reads[j])>i and reads[j][i]==N]) for N in ['A','C','T','G']})
        #m = max(c[i].values())
        m = max(c.values())
        l = [N for N in ['T', 'G', 'C', 'A'] if c[N] == m]
        bestseqs = [s + [N] for N in l for s in bestseqs]
    return "".join(bestseqs[0])


def create_fasta_per_cell_new(fastq1, fastq2, filtered_in_mapped_barcodes, output_dir,total_reads,log_file):
    with open(os.path.join(output_dir,"umi_status.csv"),'w') as log_stats:
        log_stats.write("Well_ID,cell_name,consensus_umi_barcode,read,freq,consensus_read\n")
    fastq1_dest = os.path.join(output_dir, os.path.basename(fastq1).split(".gz")[0])
    gunzip_fastq(fastq1, fastq1_dest)
    fastq2_dest = os.path.join(output_dir, os.path.basename(fastq2).split(".gz")[0])
    gunzip_fastq(fastq2, fastq2_dest)
    final_output = pd.DataFrame(columns=["Well_ID", "cell_name",'plate_total_reads', "reads_freq","umi_count", "umi_distribution"])
    with open(fastq1_dest) as f1:
        r1 = f1.readlines()
        for cell_name, cell_group in filtered_in_mapped_barcodes.groupby(by="well_coordinates"):
            well_id = cell_group["Well_ID"].iloc[0]
            cell_fasta_file = os.path.join(os.path.join(output_dir, "cells"), cell_name + ".fasta")
            with open(cell_fasta_file, 'a+') as fa:
                for consensus_umi_barcode,group in cell_group.groupby(by="consensus_umi_barcode"):
                    reads = []
                    for index,row in group.iterrows():
                        cell_barcode = row["cell_barcode"]
                        umi_barcode = row["umi_barcode"]
                        lines = subprocess.getoutput(
                        """grep -n %s %s | cut -d ":" -f1 """ % (
                        cell_barcode + umi_barcode,fastq2_dest))
                        lines = lines.split("\n")
                        lines = [int(i) - 1 for i in lines]
                        reads.extend(lines)
                    # generate a data frame of reads and their hyper variable region
                    df2 = pd.DataFrame([r1[i].rstrip() for i in reads if len(r1[i]) > 130 and "N" not in r1[i]])
                    if len(df2) == 0:
                        continue
                    df2.columns = ["read"]
                    df2_reads = df2.groupby("read").size()
                    consensus_read = df2_reads.argmax()
                    #consensus_read = consensus([list(c) for c in df2_reads.index],df2_reads.values)
                    for i, j in df2_reads.iteritems():
                        stats = well_id + "," + cell_name + "," + consensus_umi_barcode + "," + str(i[0:len(i)-1]) + "," + str(j) + "," + consensus_read
                        with open(os.path.join(output_dir, "umi_status.csv"), 'a+') as log_stats:
                            log_stats.write(stats + "\n")
                    query_line = ">" + consensus_umi_barcode + "-" + str(df2_reads.sum()) + "\n"
                    fa.write(query_line)
                    fa.write(consensus_read + "\n")
                    if well_id not in final_output.index:
                        final_output.loc[well_id,["Well_ID","cell_name",'plate_total_reads']] =  [well_id,cell_name,total_reads]
                        final_output.loc[well_id, "reads_freq"] = int(len(df2))
                        final_output.loc[well_id, "umi_count"] = 1
                        final_output.loc[well_id,"umi_distribution"] = str(int(len(df2)))
                    else:
                        final_output.loc[well_id, "reads_freq"] += int(len(df2))
                        final_output.loc[well_id, "umi_count"] += 1
                        final_output.loc[well_id, "umi_distribution"] += "," + str(int(len(df2)))
                    final_output = final_output.sort_values(by="reads_freq", ascending=False)
                    final_output.to_csv(os.path.join(output_dir, "final_output.csv"), index=False)
    with open(log_file, 'a+') as log:
        log.write(str(final_output["reads_freq"].sum()) + ",")
        log.write(str(final_output["umi_count"].sum()) + ",")
        log.write(str(final_output["Well_ID"].nunique()) + "\n")
    os.remove(fastq1_dest)
    os.remove(fastq2_dest)




def create_fasta_per_cell(fastq1, fastq2, filtered_in_mapped_barcodes, output_dir,total_reads,log):
    fastq1_dest = os.path.join(output_dir, os.path.basename(fastq1).split(".gz")[0])
    gunzip_fastq(fastq1, fastq1_dest)
    fastq2_dest = os.path.join(output_dir, os.path.basename(fastq2).split(".gz")[0])
    gunzip_fastq(fastq2, fastq2_dest)
    final_output = pd.DataFrame(columns=["Well_ID", "cell_name",'plate_total_reads', "reads_freq", "umi_distribution", "unique_var_region"])
    with open(fastq1_dest) as f1:
        r1 = f1.readlines()
        for name,group in filtered_in_mapped_barcodes.groupby(by="well_coordinates"):
            cell_name = name
            cell_fasta_file = output_dir + "/" + cell_name + ".fasta"
            reads = []
            umi_dict = dict()
            for index,row in group.iterrows():
                cell_barcode = row["cell_barcode"]
                umi_barcode = row["umi_barcode"]
                consensus_umi_barcode = row["consensus_umi_barcode"]
                lines = subprocess.getoutput(
                """grep -n %s %s | cut -d ":" -f1 """ % (
                cell_barcode + umi_barcode,fastq2_dest))
                lines = lines.split("\n")
                lines = [int(i) - 1 for i in lines]
                umi_dict.update({key:consensus_umi_barcode for key in lines})
                reads.extend(lines)
            # generate a data frame of reads and their hyper variable region
            df2 = pd.DataFrame([(r1[i], r1[i][80:130],umi_dict[i]) for i in reads if len(r1[i]) > 130])
            if len(df2) == 0:
                continue
            df2.columns = ["read","hvr","consensus_umi_barcode"]
            df_full_reads = pd.DataFrame(df2.groupby(by=["hvr","consensus_umi_barcode", "read"]).size().sort_values(ascending=False))
            df_full_reads.columns = ["read_freq"]
            df_full_reads.reset_index(inplace=True)
            df_hvr = pd.DataFrame(df2.groupby(by=["hvr"]).size().sort_values(ascending=False).head())
            df_hvr.columns = ["hvr_freq"]
            df_hvr.reset_index(inplace=True)
            # peeking only reads with abundant hyper variable region
            m = pd.merge(df_full_reads, df_hvr, on="hvr", how='right').sort_values(by="read_freq",ascending=False).reset_index(drop=True)
            m = m[~m["read"].str.contains("N")].reset_index(drop=True)
            if len(m) == 0:
                continue
            with open(cell_fasta_file, 'a') as fa:
                for index,row in m.iterrows():
                    query_line = ">" + str(index) + ":" + row["consensus_umi_barcode"] + "-" + str(row["read_freq"]) + "\n"
                    fa.write(query_line)
                    fa.write(row["read"])
            well_id = group["Well_ID"].iloc[0]
            final_output = final_output.append(
                    [{"Well_ID": well_id, "cell_name": cell_name, 'plate_total_reads': total_reads,
                      "reads_freq": m["read_freq"].sum(),
                      "umi_distribution": " ".join([str(int(count)) for count in
                                                    m.groupby("consensus_umi_barcode")["read_freq"].sum().sort_values(
                                                        ascending=False).values]),"unique_var_region": " ".join([str(int(count)) for count in
                                                    m.groupby("hvr")["read_freq"].sum().sort_values(
                                                        ascending=False).values])}])
            final_output = final_output.sort_values(by="reads_freq", ascending=False)
            final_output.to_csv(os.path.join(output_dir, "final_output.csv"), index=False)
    log.write("total_filterd_mapped_reads_after_unique_reads_filtering\t" + str(final_output["reads_freq"].sum()) + "\n")
    log.write("total_filterd_mapped_wells_after_unique_reads_filtering\t" + str(final_output["Well_ID"].nunique()) + "\n")
    os.remove(fastq1_dest)
    os.remove(fastq2_dest)


def split_to_cells(plate_name,wells_cells_file,output_dir,fastq1,fastq2,f):
    ## In order to reduce the noise and running time, we need to drop out unwanted reads
    # first way to handle it: filter by read2 (by kmers):
    ## 1. drop out reads with kmers frequency less than 4
    ## 2. droup out reads with same (similar or identitcal) umi barcode but different cell barcode and keep only the reads with the abundant cell barcode
    # second way to handle it: filter by read1 (by the gene sequene)
    ##
    log_file = os.path.join(output_dir,"split_log.csv")
    column = subprocess.getoutput(
        """gunzip -c %s | awk 'NR%s==2' | cut -b 1-15 | sort | uniq -c | sort -n """ % (fastq2, "%4")).split("\n")
    columns = [(int(column[i].strip().split(" ")[0]), column[i].strip().split(" ")[1][0:7],
                column[i].strip().split(" ")[1][7:15]) for i in range(0, len(column))]
    plate_seqs = pd.DataFrame(columns, columns=["num", "cell_barcode", "umi_barcode"])
    total_reads = plate_seqs["num"].sum()
    t = plate_seqs["num"].quantile(f)
    #t = 2
    # reading wells cells file (mapping from cell barcode to well id/well coordinates
    map_cell_to_barcode = pd.read_csv(wells_cells_file, delimiter='\t',usecols = ['Well_ID','well_coordinates', 'Cell_barcode', 'Amp_batch_ID'])
    map_cell_to_barcode = map_cell_to_barcode.loc[map_cell_to_barcode['Amp_batch_ID'] == plate_name,]
    m = pd.merge(plate_seqs, map_cell_to_barcode, left_on="cell_barcode", right_on="Cell_barcode", how='outer')
    del plate_seqs
    m = m.sort_values(by="num", ascending=False).reset_index(drop=True)
    m["consensus_umi_barcode"] = None
    mapped_barcodes = copy.deepcopy(m.loc[m["Well_ID"].notnull(),])
    un_mapped_barcodes = copy.deepcopy(m.loc[m["Well_ID"].isnull(),])

    # generating barcode to well dict
    barcode_dict, well_to_barcode = get_barcode_dict(map_cell_to_barcode,un_mapped_barcodes["cell_barcode"].unique())# generate a merged table
    del map_cell_to_barcode
    quantile = ";".join([str(x) + ":" + str(int(m["num"].quantile(x))) for x in np.append(np.linspace(0.1,0.8,8),np.linspace(0.9, 0.99, 10))])
    thresh = mapped_barcodes[mapped_barcodes["well_coordinates"].isin(["O1", "O2", "P1", "P2"])]["num"].max()
    un_mapped_barcodes = copy.deepcopy(un_mapped_barcodes.loc[un_mapped_barcodes["num"] >= t,])
    filtered_in_mapped_barcodes = copy.deepcopy(mapped_barcodes.loc[mapped_barcodes["num"] >= t,])

    with open(log_file, 'w+') as log:
        log.write("total_reads,cut_quantile,kmers_threshold,quantile,mapped_reads_before_cut,unmapped_reads_before_cut,"
                  "mapped_wells_before_cut,mapped_kmers_before_cut,unmapped_kmers_before_cut,control_threshold_before_cut,"
                  "mapped_reads_after_cut, unmapped_reads_after_cut,mapped_wells_after_cut,mapped_kmers_after_cut,unmapped_kmers_after_cut,"
                  "mapped_reads_after_unmapped_addition, unmapped_reads_after_unmapped_addition,mapped_wells_after_unmapped_addition,"
                  "mapped_kmers_after_unmapped_addition,unmapped_kmers_after_unmapped_addition,"
                  "mapped_reads_after_filtering,mapped_wells_after_filtering,mapped_kmers_after_filtering,control_threshold_after_filtering,"
                  "mapped_reads_after_wells_filtering, mapped_umis_after_wells_filtering, mapped_wells_after_wells_filtering\n")
        log.write(str(total_reads) + "," + str(f) + "," + str(t) + ",")
        log.write(quantile + ",")
        log.write(str(mapped_barcodes["num"].sum()) + "," + str(un_mapped_barcodes["num"].sum()) + ",")
        log.write(str(mapped_barcodes["Well_ID"].nunique()) + ",")
        log.write(str(len(mapped_barcodes)) + ",")
        log.write(str(len(un_mapped_barcodes)) + ",")
        if thresh > 0:
            log.write(str(thresh) +",")
        else:
            log.write(str(0) + ",")
        log.write(str(filtered_in_mapped_barcodes["num"].sum()) + ",")
        log.write(str(un_mapped_barcodes["num"].sum()) + ",")
        log.write(str(filtered_in_mapped_barcodes["Well_ID"].nunique()) + ",")
        log.write(str(len(filtered_in_mapped_barcodes)) + ",")
        log.write(str(len(un_mapped_barcodes)) + ",")

    del mapped_barcodes

    # iterating on unmapped reads with cell barcode similar to mapped one
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
    del sorted_wells
    del barcode_dict
    del well_to_barcode
    filtered_in_mapped_barcodes[["Well_ID","well_coordinates","num","cell_barcode","umi_barcode"]].to_csv(os.path.join(output_dir,"high_conf_before_filtering.csv"),index=False)
    filtered_in_mapped_barcodes = filtered_in_mapped_barcodes.sort_index()

    with open(log_file, 'a+') as log:
        log.write(str(filtered_in_mapped_barcodes["num"].sum()) + ",")
        log.write(str(un_mapped_barcodes["num"].sum()) + ",")
        log.write(str(filtered_in_mapped_barcodes["Well_ID"].nunique()) + ",")
        log.write(str(len(filtered_in_mapped_barcodes)) + ",")
        log.write(str(len(un_mapped_barcodes)) + ",")

    del un_mapped_barcodes
    # consesnus sequence for each umi barcode
    unique_consensuses = []
    for i, row in filtered_in_mapped_barcodes.iterrows():
        found = 0
        for cons_umi in unique_consensuses:
            if hamming_distance(row["umi_barcode"], cons_umi,3) <= 2:
                filtered_in_mapped_barcodes.loc[i, "consensus_umi_barcode"] = cons_umi
                found = 1
                break
        if found == 0:
            unique_consensuses.append(row["umi_barcode"])
            filtered_in_mapped_barcodes.loc[i, "consensus_umi_barcode"] = row["umi_barcode"]


    stats = filtered_in_mapped_barcodes.groupby(["consensus_umi_barcode", "Well_ID", "cell_barcode"])[
                "num"].sum().sort_index()
    stats.to_csv(os.path.join(output_dir,"molecules_stats_after_cut.csv"),header=True)

    filtered_in_mapped_barcodes[["Well_ID","well_coordinates","num","cell_barcode","umi_barcode","consensus_umi_barcode"]].to_csv(os.path.join(output_dir,"high_conf_before_filtering_after_consensus.csv"),index=False)


    # for each umi in each well - decide what is the original sequence
    """
    filtered_in_mapped_barcodes = filtered_in_mapped_barcodes.sort_index()
    for name, group in filtered_in_mapped_barcodes.groupby(by="Well_ID"):
        if group["umi_barcode"].nunique()>1:
            dom_umi = [group.iloc[0]["umi_barcode"]]
            for index,row in group.iterrows():
                for dom in dom_umi:
                    if hamming_distance(row["umi_barcode"],dom,3) <=2:
                        filtered_in_mapped_barcodes.loc[index,"consensus_umi_barcode"] = dom
                        break
                if filtered_in_mapped_barcodes.loc[index,"consensus_umi_barcode"] is None:
                    dom_umi.append(row["umi_barcode"])
                    filtered_in_mapped_barcodes.loc[index,"consensus_umi_barcode"] = row["umi_barcode"]
        else:
            for index, row in group.iterrows():
                filtered_in_mapped_barcodes.loc[index,"consensus_umi_barcode"] = row["umi_barcode"]
    filtered_in_mapped_barcodes = filtered_in_mapped_barcodes.sort_index()
    
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
    """

    # dropping out kmers with identitcal consensus umi sequence, only the most abundant well survives

    field = "consensus_umi_barcode"
    for name,group in filtered_in_mapped_barcodes.groupby(by=field):
        if group["Well_ID"].nunique()>1:
            dom_well = group.iloc[0]["Well_ID"]
            filtered_in_mapped_barcodes = filtered_in_mapped_barcodes.drop(group[group["Well_ID"] != dom_well].index)

    filtered_in_mapped_barcodes[["Well_ID","well_coordinates","num","cell_barcode","umi_barcode","consensus_umi_barcode"]].to_csv(os.path.join(output_dir,"high_conf_before_control_cut.csv"),index=False)
    old_thresh = 0 if thresh <= 0 else thresh
    thresh = filtered_in_mapped_barcodes[filtered_in_mapped_barcodes["well_coordinates"].isin(["O1", "O2", "P1", "P2"])]["num"].max()

    with open(log_file, 'a+') as log:
        log.write(str(filtered_in_mapped_barcodes["num"].sum()) + ",")
        log.write(str(filtered_in_mapped_barcodes["Well_ID"].nunique()) + ",")
        log.write(str(len(filtered_in_mapped_barcodes)) + ",")
        if thresh > 0:
            log.write(str(thresh) +",")
        else:
            log.write(str(old_thresh) + ",")
    create_fasta_per_cell_new(fastq1, fastq2, filtered_in_mapped_barcodes, output_dir,total_reads,log_file)

