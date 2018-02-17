import os
import pytest
from sctools.metrics import Runner
import tempfile
import posixpath
import shutil


data_dir = os.path.split(__file__)[0] + '/data/'
bam_file = data_dir + 'test_optimus.bam'
testing_directory = tempfile.mkdtemp()


@pytest.fixture(scope='module')
def metrics_runner():
    return Runner(bam_file, 'CR', 'UR', posixpath.join(testing_directory, 'sctools_metrics'))


def test_gather_function(metrics_runner):
    metrics_runner.run_metrics(['_noop'])


# clean up
shutil.rmtree(testing_directory)
