import urllib

import buildPlotWithDomains as domainBuilder
import dataparser as parser
import helperFiles.buildplot as plotBuilder

samples = parser.samples

for category in samples:
    for gene in samples[category]:
        chromosome, start, end = parser.getStartandEnd(gene)
        start = plotBuilder.lift(start, chromosome)
        end = plotBuilder.lift(end, chromosome)

        url = 'https://www.uniprot.org/uniprot/'

        params = {
            'query': f'gene:{gene} organism:human',
            'format': 'list',
            'columns': 'id',
            'limit': 1
        }
        data = urllib.parse.urlencode(params)
        print(data)
        data = data.encode('utf-8')
        req = urllib.request.Request(url, data)
        with urllib.request.urlopen(req) as f:
            response = f.read().decode('utf-8')
        acc = response.split()[0]  # gets the accession of the genename

        ranges = domainBuilder.buildRanges(gene, acc, chromosome)
        print(gene)
        print("hereeeeeee")
        print(ranges)
        quit()
