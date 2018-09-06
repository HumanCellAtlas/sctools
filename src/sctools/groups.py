"""
Group QC outputs

"""

from crimson import picard
import argparse
import os
import pandas as pd


def aggregate_picard_metrics_by_row(file_names, output_name):
    """
    pipeline output picard QC metrics at sinle cell/sample level.
    This functin is called to merge/aggregate QC metrics by metrics type and
    then merge multiple QC measurement into single matrix file.
    In this file, column is sample/cell and row is QC metrics
    :param files: metric files from pipeline outputs
    :param met_name: metrics name with workflow name and subworkflow name
    as prefix.such as 'run_pipelines.RunStarPipeline.alignment_summary_metrics'
    """
    # initial output
    metrics = {}
    d = pd.DataFrame()
    for file_name in file_names:
        cell_id = os.path.basename(file_name).split('_qc')[0]
        metrics[cell_id] = {}
        parsed = picard.parse(file_name)
        class_name = parsed['metrics']['class'].split('.')[2]
        # Aignment metrics return multiple lines,
        # but only output PAIRED-READS/third line
        contents = parsed['metrics']['contents']
        if class_name == "AlignmentSummaryMetrics":
            # parse out PE, R1 and R2
            rows = {}
            for m in contents:
                cat = m['CATEGORY']
                rows.update({k + '.' + cat: v for k, v in m.items() if k not in ['SAMPLE', 'LIBRARY', 'READ_GROUP', 'CATEGORY']})
        # sometimes(very rare), insertion metrics also return multiple lines
        # results to include TANDEM repeats. but we only output the first line.
        elif class_name == "InsertSizeMetrics":
            # if the elemnet counts is less than 21,
            # it means insertion metrics returns multiple line results.
            if len(contents) < 21:
                rows = contents[0]
            else:
                rows = contents
        else:
            # other metrics(so far) only return one line results.
            rows = contents
        metrics[cell_id].update({
                k: rows[k] for k in rows if k not in
                ['SAMPLE', 'LIBRARY', 'READ_GROUP', 'CATEGORY']
                })
        df = pd.DataFrame.from_dict(metrics, orient='columns')
        df.insert(0, 'Class', class_name)
        d = d.append(df)
    d.to_csv(output_name + '.csv')


def aggregate_picard_metrics_by_table(file_names, output_name):
    for file_name in file_names:
        cell_id = os.path.basename(file_name).split('_qc')[0]
        class_name = os.path.basename(file_name).split('.')[1]
        parsed = picard.parse(file_name)
        dat = pd.DataFrame.from_dict(parsed['metrics']['contents'])
        dat.insert(0, 'Sample', cell_id)
        dat.to_csv(output_name + "_" + class_name + '.csv', index=False)


def aggregate_qc_metrics(file_names, output_name):
    df = pd.DataFrame()
    for file_name in file_names:
        dat = pd.read_csv(file_name)
        df = df.append(dat, ignore_index=True)
    df.rename(columns={'Unnamed: 0': 'Metrics'}, inplace=True)
    df.to_csv(output_name + '.csv', index=False)


def parse_hisat2_log(file_names, output_name):
    metrics = {}
    for file_name in file_names:
        if '_qc' in file_name:
            cell_id = os.path.basename(file_name).split('_qc')[0]
            tag = "HISAT2G"
        elif '_rsem' in file_name:
            cell_id = os.path.basename(file_name).split('_rsem')[0]
            tag = "HISAT2T"
        with open(file_name) as f:
            dat = f.readlines()
            d = [x.strip().split(':') for x in dat]
            d.pop(0)
            metrics[cell_id] = {x[0]: x[1].strip().split(' ')[0] for x in d}
    df = pd.DataFrame.from_dict(metrics, orient='columns')
    df.insert(0, "Class", tag)
    df.to_csv(output_name + '.csv')


def parse_rsem_cnt(file_names, output_name):
    metrics = {}
    for file_name in file_names:
        cell_id = os.path.basename(file_name).split('_rsem')[0]
        i = 0
        with open(file_name) as f:
            while i < 3:
                if i == 0:
                    [N0, N1, N2, N_tot] = f.readline().strip().split(" ")
                elif i == 1:
                    [nUnique, nMulti, nUncertain] = \
                        f.readline().strip().split(" ")
                elif i == 2:
                    [nHits, read_type] = f.readline().strip().split(" ")
                i = i+1
        metrics[cell_id] = {
                "unalignable reads": N0,
                "alignable reads": N1,
                "filtered reads": N2,
                "total reads": N_tot,
                "unique aligned": nUnique,
                "multiple mapped": nMulti,
                "total alignments": nHits,
                "strand": read_type,
                "uncertain reads": nUncertain
        }
    df = pd.DataFrame.from_dict(metrics, orient='columns')
    df.insert(0, "Class", "RSEM")
    df.to_csv(output_name + '.csv')
