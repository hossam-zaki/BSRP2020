import ensembl_rest

json = ensembl_rest.symbol_lookup(
    species='homo sapiens',
    symbol='RAD51B',
    params={'expand': True}
)

start = json['start']
end = json['end']
print(end-start)


def getTranscriptsAndRanges():
    transcriptAndRange = {}
    for transcript in json['Transcript']:
        transcriptAndRange[transcript['id'] + "/" + transcript['biotype']] = (
            (transcript['start']-start, transcript['end'] - start))
    return transcriptAndRange
