import subprocess
from readtagger.cap3 import Cap3Assembly


def test_cap3_execption_handling(mocker):  # noqa: D103
    sequences = ['ATGCATGATGC']

    def raises_exception(*args, **kwargs):
        raise subprocess.CalledProcessError(127, 'Oops')

    mocker.patch('subprocess.check_call', raises_exception)
    assert len(Cap3Assembly(sequences).contigs) == 0
