from readtagger.cluster import non_evidence

NON_SUPPORT = 'non_support_test.bam'


def test_nonevidence(datadir):  # noqa: D103
    input_path = datadir[NON_SUPPORT]
    chunk = [(0, 5095, 5096, {'HWI-D00405:129:C6KNAANXX:3:1204:6411:36177',
                              'HWI-D00405:129:C6KNAANXX:3:1306:11387:59031',
                              'HWI-D00405:129:C6KNAANXX:4:2107:11019:16386',
                              'HWI-D00405:129:C6KNAANXX:3:2102:18872:29300',
                              'HWI-D00405:129:C6KNAANXX:4:1107:2350:7034',
                              'HWI-D00405:129:C6KNAANXX:4:1203:4769:64885',
                              'HWI-D00405:129:C6KNAANXX:4:2211:6930:94448'}),
             (1, 6056, 6057, {'HWI-D00405:129:C6KNAANXX:4:1309:8446:97215'}),
             (2, 16008, 16009, {'HWI-D00405:129:C6KNAANXX:4:1113:5245:46615'}),
             (3, 16422, 16423, {'HWI-D00405:129:C6KNAANXX:4:1310:9444:11374'})]
    data = {'chromosome': '2L', 'input_path': input_path, 'chunk': chunk}
    result = non_evidence(data)
    assert len(result) == 4
    assert result[0] == 12
    assert result[1] == 38
    assert result[2] == 31
    assert result[3] == 30
