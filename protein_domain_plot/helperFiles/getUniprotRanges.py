import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET


def getXML(gene):
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
    print(acc)

    url = f'https://www.uniprot.org/uniprot/{acc}.xml'
    print(url)
    req = urllib.request.Request(url)
    with urllib.request.urlopen(req) as f:
        return f.read()
        # with open(f'{acc}.xml', 'wb') as xml:
        #     xml.write(f.read())


def getRange(gene):
    root = ET.fromstring(getXML(gene))
    rangesWithLabels = {}
    for feature in root.findall("{http://uniprot.org/uniprot}entry/{http://uniprot.org/uniprot}feature"):
        beginAttr = feature.find(
            "{http://uniprot.org/uniprot}location/{http://uniprot.org/uniprot}begin")
        if beginAttr is not None:
            begin = beginAttr.attrib['position']
            end = feature.find(
                "{http://uniprot.org/uniprot}location/{http://uniprot.org/uniprot}end").attrib['position']
            if(feature.attrib['type'] != 'domain' and 'region' not in feature.attrib['type']):
                continue
            description = feature.attrib['description']
            rangesWithLabels[description] = ((begin, end))
    return rangesWithLabels


if __name__ == "__main__":
    print(getRange('XRCC6'))
