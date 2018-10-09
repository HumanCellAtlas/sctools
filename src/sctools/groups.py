"""
Group QC outputs

"""

from crimson import picard
import os
import pandas as pd


def write_aggregated_picard_metrics_by_row(file_names, output_name):
    """Command line entrypoint to parse, aggreagete and write Picard row metrics.
    Parameters
    ----------
    args:
        file_names: array of files. the basename of inputs should be formated
        as 'samplename_qc',such as
        "samplename_qc.alignment_summary_metrics.txt" and "samplename_qc.insert_size_metrics.txt"
        output_name: prefix of output file name without extension.
    Returns
    ----------
        return: 0
        return if the program completes successfully.
    """
    # initial output
    metrics = {}
    d = pd.DataFrame()
    for file_name in file_names:
        cell_id = os.path.basename(file_name).split('_qc')[0]
        metrics[cell_id] = {}
        parsed = picard.parse(file_name)
        class_name = parsed['metrics']['class'].split('.')[2]
        # Alignment metrics return multiple lines,
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
            # if the element counts is less than 21,
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
    d_T = d.T
    d_T.to_csv(output_name + '.csv')


def write_aggregated_picard_metrics_by_table(file_names, output_name):
    """Command line entrypoint to parse and write Picard table metrics.
    Parameters
    ----------
    args:
        file_names: array of files.the basename of inputs should be formated as 'samplename_qc'
        output_name: prefix of output file name. the basename of outputs
        includes the Picard metrics class name.
    Returns
    ----------
        return: 0
        return if the program completes successfully.
    """
    for file_name in file_names:
        cell_id = os.path.basename(file_name).split('_qc')[0]
        class_name = os.path.basename(file_name).split('.')[1]
        parsed = picard.parse(file_name)
        dat = pd.DataFrame.from_dict(parsed['metrics']['contents'])
        dat.insert(0, 'Sample', cell_id)
        dat.to_csv(output_name + "_" + class_name + '.csv', index=False)


def write_aggregated_qc_metrics(file_names, output_name):
    """Command line entrypoint to merge Picard metrics along with RSEM and HISAT2 log
    Parameters
    ----------
    args:
        file_names: array of files,such as Picard row metric, hisat2 metrics.
        output_name: prefix of output file name.
    Returns
    ----------
        return: 0
        return if the program completes successfully.
    """
    df = pd.DataFrame()
    for file_name in file_names:
        dat = pd.read_csv(file_name, index_col=0)
        print(dat.index)
        df = pd.concat([df, dat], axis=1, join_axes=[dat.index])
    df.to_csv(output_name + '.csv', index=True)


def parse_hisat2_log(file_names, output_name):
    """Command line entrypoint parse, aggreagete and write HISAT2 logs
    Parameters
    ----------
    args:
        file_names: array of HISAT2 log files. Basename of file indicates the alignment references
        '_qc' indicates the genome reference and '_rsem' indicates the transcriptome reference alignment.
        output_name: prefix of output file name. 
    Returns
    ----------
        return: 0
        return if the program completes successfully.
    """
    metrics = {}
    tag = "NONE"
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
            # remove the first row of each section.
            d.pop(0)
            metrics[cell_id] = {x[0]: x[1].strip().split(' ')[0] for x in d}
    df = pd.DataFrame.from_dict(metrics, orient='columns')
    df.insert(0, "Class", tag)
    df_T = df.T
    df_T.to_csv(output_name + '.csv')


def parse_rsem_cnt(file_names, output_name):
    """Command line entrypoint parse, aggreagete and write RSEM cnt
    Parameters
    ----------
    args:
        file_names: array of RSEM cnt files. The basename of inputs should be '_rsem'
        output_name: prefix of output file name.
    Returns
    ----------
        return: 0
        return if the program completes successfully.
    """
    metrics = {}
    for file_name in file_names:
        cell_id = os.path.basename(file_name).split('_rsem')[0]
        i = 0
        with open(file_name) as f:
            while i < 3:
                if i == 0:
                    [N0, N1, N2, N_tot] = f.readline().strip().split(" ")
                elif i == 1:
                    [n_unique, n_multi, n_uncertain] = \
                        f.readline().strip().split(" ")
                elif i == 2:
                    [n_hits, read_type] = f.readline().strip().split(" ")
                i = i+1
        metrics[cell_id] = {
                "unalignable reads": N0,
                "alignable reads": N1,
                "filtered reads": N2,
                "total reads": N_tot,
                "unique aligned": n_unique,
                "multiple mapped": n_multi,
                "total alignments": n_hits,
                "strand": read_type,
                "uncertain reads": n_uncertain
        }
    df = pd.DataFrame.from_dict(metrics, orient='columns')
    df.insert(0, "Class", "RSEM")
    df_T = df.T
    df_T.to_csv(output_name + '.csv')
