import ensembl_rest

json = ensembl_rest.symbol_lookup(
    species='homo sapiens',
    symbol='RAD51B',
    params={'expand': True}
)


start = json['start']
end = json['end']


def getTranscriptsAndRanges():
    transcriptAndRange = {}
    for transcript in json['Transcript']:
        if transcript['biotype'] != 'protein_coding':
            continue
        transcriptAndRange[transcript['id']] = []
        for exon in transcript['Exon']:
            transcriptAndRange[transcript
                               ['id']].append((exon['start']-start, exon['end']-start))
    print(transcriptAndRange)
    return transcriptAndRange
