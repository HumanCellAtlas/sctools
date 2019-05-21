import os
import csv
import itertools
from collections import OrderedDict
from sctools import platform


data_dir = os.path.split(__file__)[0] + '/data/group_metrics/'
unpaired_data_dir = os.path.split(__file__)[0] + '/data/group_metrics_unpaired_ss2/'


def check_parsed_metrics_csv(file_name, cell_id, class_name, expected_metrics):
    with open(file_name) as f:
        column_headers = f.readline().strip().split(',')
        classes = f.readline().strip().split(',')
        metrics = f.readline().strip().split(',')
        assert classes[0] == 'Class'
        assert set(classes[1:]) == {class_name}
        for idx, each in enumerate(column_headers):
            if idx == 0:
                assert metrics[0] == cell_id
            if idx > 0:
                metric_name = column_headers[idx]
                assert metrics[idx] == expected_metrics[metric_name]


def test_write_aggregated_picard_metrics_by_row():
    args = [
        '-f',
        data_dir + 'test_qc.alignment_summary_metrics.txt',
        data_dir + 'test_qc.insert_size_metrics.txt',
        data_dir + 'test_qc.duplicate_metrics.txt',
        data_dir + 'test_qc.rna_metrics.txt',
        data_dir + 'test_qc.gc_bias.summary_metrics.txt',
        '-t',
        'Picard',
        '-o',
        'output_picard_group',
    ]
    return_code = platform.GenericPlatform.group_qc_outputs(args)
    assert return_code == 0

    expected_metrics = {}
    with open(data_dir + 'expected_picard_group.csv') as f:
        column_headers = f.readline().strip().split(',')
        classes = f.readline().strip().split(',')
        metrics = f.readline().strip().split(',')
        for idx, each in enumerate(column_headers):
            expected_metrics[each] = {'class': classes[idx], 'metric': metrics[idx]}
    with open('output_picard_group.csv') as f:
        column_headers = f.readline().strip().split(',')
        classes = f.readline().strip().split(',')
        metrics = f.readline().strip().split(',')
        assert len(column_headers) == len(expected_metrics.keys())
        for idx, each in enumerate(column_headers):
            header = expected_metrics[each]
            assert classes[idx] == header['class']
            assert metrics[idx] == header['metric']
    os.remove('output_picard_group.csv')


def test_write_aggregated_picard_metrics_by_table():
    args = [
        '-t',
        'PicardTable',
        '-o',
        'output_picard_group',
        '-f',
        data_dir + 'test_qc.error_summary_metrics.txt',
    ]
    return_code = platform.GenericPlatform.group_qc_outputs(args)
    assert return_code == 0

    expected_metrics = [
        OrderedDict(
            [
                ('Sample', 'test'),
                ('ALT_BASE', 'C'),
                ('ALT_COUNT', '16'),
                ('REF_BASE', 'A'),
                ('REF_COUNT', '231512'),
                ('SUBSTITUTION', 'A>C'),
                ('SUBSTITUTION_RATE', '6.9e-05'),
            ]
        ),
        OrderedDict(
            [
                ('Sample', 'test'),
                ('ALT_BASE', 'G'),
                ('ALT_COUNT', '156'),
                ('REF_BASE', 'A'),
                ('REF_COUNT', '231512'),
                ('SUBSTITUTION', 'A>G'),
                ('SUBSTITUTION_RATE', '0.000673'),
            ]
        ),
        OrderedDict(
            [
                ('Sample', 'test'),
                ('ALT_BASE', 'T'),
                ('ALT_COUNT', '16'),
                ('REF_BASE', 'A'),
                ('REF_COUNT', '231512'),
                ('SUBSTITUTION', 'A>T'),
                ('SUBSTITUTION_RATE', '6.9e-05'),
            ]
        ),
        OrderedDict(
            [
                ('Sample', 'test'),
                ('ALT_BASE', 'A'),
                ('ALT_COUNT', '16'),
                ('REF_BASE', 'C'),
                ('REF_COUNT', '173880'),
                ('SUBSTITUTION', 'C>A'),
                ('SUBSTITUTION_RATE', '9.2e-05'),
            ]
        ),
        OrderedDict(
            [
                ('Sample', 'test'),
                ('ALT_BASE', 'G'),
                ('ALT_COUNT', '14'),
                ('REF_BASE', 'C'),
                ('REF_COUNT', '173880'),
                ('SUBSTITUTION', 'C>G'),
                ('SUBSTITUTION_RATE', '8.1e-05'),
            ]
        ),
        OrderedDict(
            [
                ('Sample', 'test'),
                ('ALT_BASE', 'T'),
                ('ALT_COUNT', '82'),
                ('REF_BASE', 'C'),
                ('REF_COUNT', '173880'),
                ('SUBSTITUTION', 'C>T'),
                ('SUBSTITUTION_RATE', '0.000471'),
            ]
        ),
    ]

    with open('output_picard_group_error_summary_metrics.csv') as f:
        reader = csv.DictReader(f)
        for line in reader:
            assert line in expected_metrics
    os.remove('output_picard_group_error_summary_metrics.csv')


def test_parse_hisat2_paired_end_log():
    args = [
        '-f',
        data_dir + 'test_hisat2_paired_end_qc.log',
        '-t',
        'HISAT2',
        '-o',
        'output_hisat2',
    ]
    return_code = platform.GenericPlatform.group_qc_outputs(args)
    assert return_code == 0

    cell_id = 'test_hisat2_paired_end'
    tag = 'HISAT2G'
    expected_metrics = {
        'Total pairs': '5479',
        'Aligned concordantly or discordantly 0 time': '412',
        'Aligned concordantly 1 time': '4414',
        'Aligned concordantly >1 times': '652',
        'Aligned discordantly 1 time': '1',
        'Total unpaired reads': '824',
        'Aligned 0 time': '478',
        'Aligned 1 time': '240',
        'Aligned >1 times': '106',
        'Overall alignment rate': '95.64%',
    }
    check_parsed_metrics_csv('output_hisat2.csv', cell_id, tag, expected_metrics)
    os.remove('output_hisat2.csv')


def test_parse_hisat2_transcriptome_log():
    args = [
        '-f',
        data_dir + 'test_hisat2_transcriptome_rsem.log',
        '-t',
        'HISAT2',
        '-o',
        'output_hisat2_trans',
    ]
    return_code = platform.GenericPlatform.group_qc_outputs(args)
    assert return_code == 0

    cell_id = 'test_hisat2_transcriptome'
    tag = 'HISAT2T'
    expected_metrics = {
        'Total pairs': '5479',
        'Aligned concordantly or discordantly 0 time': '3635',
        'Aligned concordantly 1 time': '360',
        'Aligned concordantly >1 times': '1484',
        'Aligned discordantly 1 time': '0',
        'Total unpaired reads': '7270',
        'Aligned 0 time': '7270',
        'Aligned 1 time': '0',
        'Aligned >1 times': '0',
        'Overall alignment rate': '33.66%',
    }
    check_parsed_metrics_csv('output_hisat2_trans.csv', cell_id, tag, expected_metrics)
    os.remove('output_hisat2_trans.csv')


def test_parse_rsem_cnt():
    file_name = data_dir + 'test_rsem.cnt'
    args = ['-f', file_name, '-t', 'RSEM', '-o', 'output_rsem']
    return_code = platform.GenericPlatform.group_qc_outputs(args)
    assert return_code == 0

    cell_id = 'test'
    class_name = 'RSEM'
    expected_metrics = None
    with open(file_name) as f:
        N0, N1, N2, N_tot = f.readline().strip().split(" ")
        n_unique, n_multi, n_uncertain = f.readline().strip().split(" ")
        n_hits, read_type = f.readline().strip().split(" ")
        expected_metrics = {
            'unalignable reads': N0,
            'alignable reads': N1,
            'filtered reads': N2,
            'total reads': N_tot,
            'unique aligned': n_unique,
            'multiple mapped': n_multi,
            'total alignments': n_hits,
            'strand': read_type,
            'uncertain reads': n_uncertain,
        }
    check_parsed_metrics_csv('output_rsem.csv', cell_id, class_name, expected_metrics)
    os.remove('output_rsem.csv')


def test_write_aggregated_qc_metrics():
    input_files = [
        data_dir + 'test_picard_group.csv',
        data_dir + 'test_hisat2.csv',
        data_dir + 'test_hisat2_trans.csv',
        data_dir + 'test_rsem.csv',
    ]
    args = [
        '-f',
        data_dir + 'test_picard_group.csv',
        data_dir + 'test_hisat2.csv',
        data_dir + 'test_hisat2_trans.csv',
        data_dir + 'test_rsem.csv',
        '-t',
        'Core',
        '-o',
        'output_QCs',
    ]
    return_code = platform.GenericPlatform.group_qc_outputs(args)
    assert return_code == 0

    expected_metrics = []
    expected_headers = []
    for input_file in input_files:
        with open(input_file) as f:
            reader = csv.DictReader(f)
            expected_headers.extend(reader.fieldnames[1:])
            for idx, line in enumerate(reader):
                if len(expected_metrics) < idx + 1:
                    expected_metrics.append(line)
                else:
                    expected_metrics[idx].update(line)
    output_headers = []
    with open('output_QCs.csv') as output_file:
        reader = csv.DictReader(output_file)
        output_headers.extend(reader.fieldnames)
        for line in reader:
            assert line in expected_metrics
    # The output file should contain all of the column headers from the input files plus the "joined column" containing row headers
    assert len(output_headers) == len(expected_headers) + 1
    os.remove('output_QCs.csv')


def test_unpaired_ss2_write_aggregated_picard_metrics_by_row():

    sources = [
        unpaired_data_dir + 'SRR6258488_qc.alignment_summary_metrics.txt',
        unpaired_data_dir + 'SRR6258488_qc.duplicate_metrics.txt',
        unpaired_data_dir + 'SRR6258488_qc.gc_bias.summary_metrics.txt',
        unpaired_data_dir + 'SRR6258488_qc.rna_metrics.txt',
    ]

    args = ['-f', *sources, '-t', 'Picard', '-o', 'output_picard_group_unpaired']
    return_code = platform.GenericPlatform.group_qc_outputs(args)
    assert return_code == 0

    expected_metrics = {}

    for source in sources:
        with open(source) as f:
            for line in f:
                if line.startswith("## METRICS CLASS"):
                    class_ = line.strip().split('\t')[1].split('.')[-1]
                    break
            labels = f.readline().strip().split("\t")
            values = f.readline().strip().split("\t")

            for label, value in itertools.zip_longest(labels, values, fillvalue=""):
                if label in ("LIBRARY", "SAMPLE", "READ_GROUP", "CATEGORY"):
                    continue
                if class_ == "AlignmentSummaryMetrics":
                    label += ".UNPAIRED"
                try:
                    value = str(float(value))
                except ValueError:
                    pass
                expected_metrics[(class_, label)] = value
    expected_metrics[("Class", "")] = "SRR6258488"

    with open('output_picard_group_unpaired.csv') as f:
        labels = f.readline().strip().split(',')
        classes = f.readline().strip().split(',')
        values = f.readline().strip().split(',')
        assert len(labels) == len(expected_metrics)

        for class_, label in expected_metrics:
            if class_ not in classes or label not in labels:
                print("!", class_, label)

        for class_, label, value in zip(classes, labels, values):
            assert (class_, label) in expected_metrics
            try:
                value = str(float(value))
            except ValueError:
                value = value
            try:
                expected_value = str(float(expected_metrics[(class_, label)]))
            except ValueError:
                expected_value = expected_metrics[(class_, label)]
            assert value == expected_value
    os.remove('output_picard_group_unpaired.csv')
