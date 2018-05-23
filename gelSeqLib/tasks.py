import argparse
import os
import pickle
import subprocess
import sys
import warnings

import pandas as pd

sys.path.insert(0, '/home/labs/amit/diklag/PycharmProjects/VDJ_Dikla/python_lsf_wrapper/')
from python_lsf_wrapper.LSF import LSF, wait_for_jobs

from configparser import ConfigParser

from gelSeqLib import plate_to_cells,VDJ_func
from gelSeqLib import io_func
from gelSeqLib import base_dir



class Task:
    base_parser = argparse.ArgumentParser(add_help=False,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    base_parser.add_argument('--ncores', '-p', metavar="<CORES>",
                             help='number of processor cores to use', type=int,
                             default=1)
    base_parser.add_argument('--config_file', '-c', metavar="<CONFIG_FILE>",
                             help='config file to use',
                             default='~/.gelseqrc')

    config = None

    def run(self):
        pass

    def get_binary(self, name):
        tool_key = name.lower() + '_path'
        user_path = None
        if self.config.has_option('tool_locations', tool_key):
            user_path = self.resolve_relative_path(
                self.config.get('tool_locations', tool_key))
        return io_func.check_binary(name, user_path)

    def read_config(self, config_file):
        # Read config file
        if not config_file:
            config_file = '~/.gelseqrc'
        config_file = os.path.expanduser(config_file)
        if not os.path.isfile(config_file):
            print(
                "Config file not found at ~/.gelseqrc. Using default gelSeq.conf in repo...")
            config_file = os.path.join(base_dir, 'gelSeq.conf')
        VDJ_func.check_config_file(config_file)
        config = ConfigParser()
        config.read(config_file)
        return config

    def resolve_relative_path(self, path):
        if not path.startswith("/"):
            base_directory = os.path.abspath(os.path.dirname(__file__))
            full_path = os.path.normpath(
                "/{}/../{}".format(base_directory, path))
        else:
            full_path = path
        return full_path

    def get_index_location(self, name):
        location = os.path.join(base_dir, 'resources', self.species, name)
        return location


    def get_resources_root(self, species):
        resources_dir = os.path.join(base_dir, 'resources')
        resources_root = os.path.join(resources_dir, species)
        return resources_root


    def get_available_species(self):
        resources_dir = os.path.join(base_dir, 'resources')
        species_dirs = next(os.walk(resources_dir))[1]
        return species_dirs


class Plate_Task(Task):

    def __init__(self, **kwargs):
        if not kwargs:
            # get list of all available species in resources

            self.parser = argparse.ArgumentParser(add_help=True,description="Process fastq files from single plate, split the reads by cell basrcodes",
                parents=[self.base_parser],
                formatter_class=argparse.ArgumentDefaultsHelpFormatter)

            self.parser.add_argument('--resume_with_existing_files', '-r',
                                help='look for existing intermediate files and use those instead of starting from scratch',
                                action="store_true")

            self.parser.add_argument('--species', '-s',
                                help='Species to use for reconstruction',
                                choices=self.get_available_species(),
                                default='Hsap')

            self.parser.add_argument('fastq1', metavar="<FASTQ1>",
                                help='first fastq file - read1')
            self.parser.add_argument('fastq2', metavar="<FASTQ2>",
                                help='second fastq file - read2')
            self.parser.add_argument('plate_name', metavar="<PLATE_NAME>",help='name of plate for file labels')
            self.parser.add_argument('output_dir', metavar="<OUTPUT_DIR>",
                                help='directory for output as <output_dir>/<plate_name>')
            self.parser.add_argument('--receptor_name',
                                help="Name of receptor to reconstruct",
                                default='TCR')
            self.parser.add_argument('--loci',
                                help="Space-separated list of loci to reconstruct for receptor",
                                default=['A', 'B'], nargs='*')
            self.parser.add_argument('--full',
                                help="Continue the full process - after splitting to cells, create new job for each cell",
                                action="store_true")
            self.parser.add_argument('--filter',
                                help="umis with more than filter reads (with respect to quantile) will be saved",type=float, default=0.98)


            args = self.parser.parse_args(sys.argv[2:])
            self.plate_name = args.plate_name
            self.fastq1 = args.fastq1
            self.fastq2 = args.fastq2
            self.ncores = str(args.ncores)
            self.species = args.species
            self.resume_with_existing_files = args.resume_with_existing_files
            self.output_dir = args.output_dir
            self.receptor_name = args.receptor_name
            self.loci = args.loci
            self.full = args.full
            self.filter = args.filter
            config_file = args.config_file


        else:
            self.plate_name = kwargs.get('plate')
            self.fastq1 = kwargs.get('fastq1')
            self.fastq2 = kwargs.get('fastq2')
            self.ncores = kwargs.get('ncores')
            self.species = kwargs.get('species')
            self.resume_with_existing_files = kwargs.get(
                'resume_with_existing_files')
            self.output_dir = kwargs.get('output_dir')
            self.receptor_name = kwargs.get('receptor_name')
            self.loci = kwargs.get('loci')
            self.full = kwargs.get('full')
            self.filter = kwargs.get('filter')
            config_file = kwargs.get('config_file')

        self.config = self.read_config(config_file)

    def run(self, **kwargs):
        # Set-up output directories
        root_output_dir = os.path.abspath(self.output_dir)
        io_func.makeOutputDir(root_output_dir)
        self.output_dir = os.path.join(root_output_dir , self.plate_name)
        io_func.makeOutputDir(self.output_dir)


        self.split_to_cells()
        #cells_dir = os.path.join(self.output_dir,"cells")
        #umis_dir = os.path.join(self.output_dir, "cells_umis")
        mapper = []
        #for dir in [cells_dir,umis_dir]:
        mapper += [_CELLrun(fasta.replace('.fasta', ''), self.output_dir + "/" + fasta, self.output_dir)
                  for fasta in os.listdir(self.output_dir) if ".fasta" in fasta]
        for job in mapper:
            job.submit_command(cpu_cores=5, memory=600, queue="new-short")
        wait_for_jobs(mapper)

        #for dir in [cells_dir,umis_dir]:
        df1 = pd.read_csv(self.output_dir + "/final_output.csv")
        list_csv = [self.output_dir + "/" + file for file in os.listdir(self.output_dir) if ".csv" in file and "output" not in file and "high" not in file]
        df2 = pd.DataFrame(columns=["cell_name","V_first","V_first_counts","V_first_avg_e_value",
                                                                          "V_second","V_second_counts","V_second_avg_e_value",
                                                                          "D_first","D_first_counts","D_first_avg_e_value",
                                                                          "D_second","D_second_counts","D_second_avg_e_value",
                                                                          "J_first", "J_first_counts","J_first_avg_e_value",
                                                                          "J_second", "J_second_counts","J_second_avg_e_value",
                                                                          "CDR3_first", "CDR3_first_counts","CDR3_first_identity","CDR3_first_unique_reads",
                                                                          "CDR3_translation_first","CDR3_translation_first_counts","CDR3_translation_first_identity","CDR_translation_first_unique_reads",
                                                                          "CDR3_second", "CDR3_second_counts","CDR3_second_identity","CDR3_second_unique_reads",
                                                                          "CDR3_translation_second","CDR3_translation_second_counts","CDR3_translation_second_identity",
                                                                          "CDR3_translation_second_unique_reads"])
        for file in list_csv:
            df2 = df2.append(pd.read_csv(file), ignore_index=True)
        #for file in list_csv:
            #os.remove(file)
        df1 = pd.merge(df1,df2,on="cell_name",how="inner",left_index=False,right_index=False)
        df1 = df1[['Well_ID','cell_name','reads_freq','plate_total_reads','umi_distribution',"unique_var_region","V_first","V_first_counts","V_first_avg_e_value",
                                                                          "V_second","V_second_counts","V_second_avg_e_value",
                                                                          "D_first","D_first_counts","D_first_avg_e_value",
                                                                          "D_second","D_second_counts","D_second_avg_e_value",
                                                                          "J_first", "J_first_counts","J_first_avg_e_value",
                                                                          "J_second", "J_second_counts","J_second_avg_e_value",
                                                                          "CDR3_first", "CDR3_first_counts","CDR3_first_identity","CDR3_first_unique_reads",
                                                                          "CDR3_translation_first","CDR3_translation_first_counts","CDR3_translation_first_identity","CDR_translation_first_unique_reads",
                                                                          "CDR3_second", "CDR3_second_counts","CDR3_second_identity","CDR3_second_unique_reads",
                                                                          "CDR3_translation_second","CDR3_translation_second_counts","CDR3_translation_second_identity",
                                                                          "CDR3_translation_second_unique_reads"]]
        df1.to_csv(self.output_dir + "/final_table.csv")


    '''
    for file in os.listdir(self.output_dir):
        if ".fasta" in file or ".err" in file or ".log" in file or ".pkl" in file:
            os.remove(self.output_dir + "/" + file)
    '''


    def split_to_cells(self):
        wells_cells_file = self.config.get('project directories','wells_cells_file')
        plate_to_cells.split_to_cells(self.plate_name, wells_cells_file, self.output_dir,
                                      self.fastq1, self.fastq2,self.filter)


class Cell_Task(Task):
    def __init__(self, **kwargs):
        if not kwargs:
            self.parser = argparse.ArgumentParser(add_help=True,
                                                  description="Reconstruct TCR sequences from RNAseq reads for a single cell",
                                                  parents=[self.base_parser],
                                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter)

            self.parser.add_argument('--resume_with_existing_files', '-r',
                                     help='look for existing intermediate files and use those instead of starting from scratch',
                                     action="store_true")

            self.parser.add_argument('--species', '-s',
                                     help='Species to use for reconstruction',
                                     choices=self.get_available_species(),
                                     default='Hsap')

            self.parser.add_argument('fasta', metavar="<FASTA>",
                                     help='fasta file')
            self.parser.add_argument('cell_name', metavar="<CELL_NAME>", help='name of cell for file labels')
            self.parser.add_argument('output_dir', metavar="<OUTPUT_DIR>",
                                     help='directory for output as <output_dir>/<cell_name>')

            self.parser.add_argument('--receptor_name',
                                         help="Name of receptor to reconstruct",
                                         default='TCR')
            self.parser.add_argument('--loci',
                                         help="Space-separated list of loci to reconstruct for receptor",
                                         default=['A', 'B'], nargs='+')

            args = self.parser.parse_args(sys.argv[2:])
            self.cell_name = args.cell_name
            self.fasta = args.fasta
            self.ncores = str(args.ncores)
            self.species = args.species
            self.resume_with_existing_files = args.resume_with_existing_files
            self.output_dir = args.output_dir
            self.receptor_name = args.receptor_name
            self.loci = args.loci
            config_file = args.config_file

        else:
            self.cell_name = kwargs.get('cell')
            self.fasta = kwargs.get('fasta')
            self.ncores = kwargs.get('ncores')
            self.species = kwargs.get('species')
            self.resume_with_existing_files = kwargs.get(
                'resume_with_existing_files')
            self.output_dir = kwargs.get('output_dir')
            self.receptor_name = kwargs.get('receptor_name')
            self.loci = kwargs.get('loci')

            config_file = kwargs.get('config_file')

        self.config = self.read_config(config_file)
        # self.locus_names = ["TCRA", "TCRB"]

    def run(self, **kwargs):
        #self.collapse_uniq_seqs()
        self.cell_name=self.cell_name.split("_")[0]
        cell = self.ig_blast()

        summary = cell.choose_recombinants()
        for locus in self.loci:
            with open(
                    "{output_dir}/{cell_name}_{receptor}_{locus}.pkl".format(
                        output_dir=self.output_dir,
                        cell_name=cell.name,
                        receptor=self.receptor_name,locus=locus), 'wb') as pf:
                pickle.dump(cell, pf, protocol=0)

        if len(summary[self.receptor_name]) != 0:
            for locus in self.loci:
                summary[self.receptor_name][locus].to_csv("{output_dir}/{cell_name}_{receptor}_{locus}.csv".format(
                            output_dir=self.output_dir,
                            cell_name=cell.name,
                            receptor=self.receptor_name, locus=locus))


    def collapse_uniq_seqs(self):
        for locus in self.loci:
            new_fasta = os.path.join(self.output_dir,self.cell_name) + "_"+ self.receptor_name + "_" + locus + "_collapsed.fasta"
            subprocess.getoutput("""fastx_collapser -i %s -o %s""" % (self.fasta,new_fasta))
            #temp_uniq = subprocess.getoutput("""cat %s | awk 'NR%s==0' | sort | uniq -c | sort -n""" % (self.fasta,"%2")).split("\n")
            """
            with open(new_fasta,'w') as fasta_file:
                for row in temp_uniq:
                    seq = row.strip().split(" ")[1].strip("\n")
                    num = row.strip().split(" ")[0].strip("\n")
                    query_name = ">" + num
                    fasta_file.write(query_name + "\n")
                    fasta_file.write(seq + "\n")
            """
            self.fasta = new_fasta



    def ig_blast(self):
        igblastn = self.get_binary('igblastn')

        # Reference data locations
        igblast_index_location = self.get_index_location('igblast_dbs')
        imgt_seq_location = self.get_index_location('raw_seqs')
        aux_file_location = self.get_index_location('gl.aux')
        igblast_seqtype = self.config.get('IgBlast_options', 'igblast_seqtype')

        # IgBlast of assembled contigs
        VDJ_func.run_IgBlast(igblastn, self.fasta, self.receptor_name, self.loci,
                             self.output_dir, self.cell_name, self.species,
                             igblast_index_location,
                             igblast_seqtype, aux_file_location, self.resume_with_existing_files)
        print()

        with warnings.catch_warnings():
            cell = io_func.parse_IgBLAST(self.fasta,self.receptor_name, self.loci,
                                         self.output_dir, self.cell_name, imgt_seq_location,
                                         self.species, 'imgt',
                                         50)
        return cell



class _CELLrun(LSF):
    def __init__(self, name, fasta, output_dir,loci='B',receptor_name='TCR', species= 'Hsap'):
        """\
        name        - cell name
        fasta       - fasta path of the cell's reads
        output_dir  - dir of plate of origin
        species     - cell species (Human=Hsap, Mouse = Mmus)

        """

        # build the alignment comand
        cell_cmd = "python3.5 " + os.path.join(base_dir,"gelSeq.py") + " cell -s "+ species +" --loci=" + loci + " --receptor_name=" + receptor_name + " " + fasta + " " + name + " " + output_dir

        self.cmd = cell_cmd

        # add error and log files for the lsf output
        self.cmd = " -o " + fasta.replace(".fasta", ".log") + " " + self.cmd
        self.cmd = " -e "+ fasta.replace(".fasta", ".err") + " " + self.cmd
