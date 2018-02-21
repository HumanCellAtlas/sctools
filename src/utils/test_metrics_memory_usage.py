import os
from memory_profiler import profile
from sctools.metrics import Runner
import tempfile
import posixpath
import shutil


data_dir = os.path.split(__file__)[0] + '/data/'
bam_file = data_dir + 'test_optimus.bam'
testing_directory = tempfile.mkdtemp()


def metrics_runner():
    return Runner(bam_file, 'CR', 'UR', posixpath.join(testing_directory, 'sctools_metrics'))


@profile
def test_gather_function():
    metrics_runner().run_metrics(['_noop'])


# clean up
shutil.rmtree(testing_directory)

if __name__ == "__main__":
    test_gather_function()
