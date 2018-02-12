from readtagger.cluster import non_evidence

NON_SUPPORT = 'non_support_test.bam'


def test_nonevidence(datadir_copy):  # noqa: D103
    input_path = str(datadir_copy[NON_SUPPORT])
    for_evidence_kwargs = {0: {""}}
    chunk = [(0, 5095, 5096, {'HWI-D00405:129:C6KNAANXX:3:1204:6411:36177',
                              'HWI-D00405:129:C6KNAANXX:3:1306:11387:59031',
                              'HWI-D00405:129:C6KNAANXX:4:2107:11019:16386',
                              'HWI-D00405:129:C6KNAANXX:3:2102:18872:29300',
                              'HWI-D00405:129:C6KNAANXX:4:1107:2350:7034',
                              'HWI-D00405:129:C6KNAANXX:4:1203:4769:64885',
                              'HWI-D00405:129:C6KNAANXX:4:2211:6930:94448'},
              for_evidence_kwargs),
             (1, 6056, 6057, {'HWI-D00405:129:C6KNAANXX:4:1309:8446:97215'},
              for_evidence_kwargs),
             (2, 16008, 16009, {'HWI-D00405:129:C6KNAANXX:4:1113:5245:46615'},
              for_evidence_kwargs),
             (3, 16422, 16423, {'HWI-D00405:129:C6KNAANXX:4:1310:9444:11374'},
              for_evidence_kwargs)]
    data = {'chromosome': '2L', 'input_path': input_path, 'chunk': chunk}
    result = non_evidence(data)['against']
    assert len(result) == 4
    assert len(result[0]) == 12
    assert len(result[1]) == 38
    assert len(result[2]) == 31
    assert len(result[3]) == 30


def test_nonevidence_exception_handling(datadir_copy):  # noqa: D103
    input_path = 'some_path_that_does_not_exists'
    for_evidence_kwargs = {0: {""}}
    chunk = [(0, 5095, 5096, {'HWI-D00405:129:C6KNAANXX:3:1204:6411:36177',
                              'HWI-D00405:129:C6KNAANXX:3:1306:11387:59031',
                              'HWI-D00405:129:C6KNAANXX:4:2107:11019:16386',
                              'HWI-D00405:129:C6KNAANXX:3:2102:18872:29300',
                              'HWI-D00405:129:C6KNAANXX:4:1107:2350:7034',
                              'HWI-D00405:129:C6KNAANXX:4:1203:4769:64885',
                              'HWI-D00405:129:C6KNAANXX:4:2211:6930:94448'},
              for_evidence_kwargs),
             (1, 6056, 6057, {'HWI-D00405:129:C6KNAANXX:4:1309:8446:97215'},
              for_evidence_kwargs),
             (2, 16008, 16009, {'HWI-D00405:129:C6KNAANXX:4:1113:5245:46615'},
              for_evidence_kwargs),
             (3, 16422, 16423, {'HWI-D00405:129:C6KNAANXX:4:1310:9444:11374'},
              for_evidence_kwargs)]
    data = {'chromosome': '2L', 'input_path': input_path, 'chunk': chunk}
    assert len(non_evidence(data)['against']) == 0
