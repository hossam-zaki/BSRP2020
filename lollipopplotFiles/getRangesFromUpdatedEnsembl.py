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
        if 'Translation' not in transcript:
            continue
        transcriptAndRange[transcript['Translation']['id']] = []
        for exon in transcript['Exon']:
            transcriptAndRange[transcript['Translation']
                               ['id']].append((exon['start']-start, exon['end']-start))
    return transcriptAndRange
