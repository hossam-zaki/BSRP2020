import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET


def getXML():
    url = 'https://www.uniprot.org/uniprot/'

    params = {
        'query': 'gene:rad51b organism:human',
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
    req = urllib.request.Request(url)
    with urllib.request.urlopen(req) as f:
        return f.read()
        # with open(f'{acc}.xml', 'wb') as xml:
        #     xml.write(f.read())


def getRange():
    root = ET.fromstring(getXML())
    rangesWithLabels = {}
    for feature in root.findall("{http://uniprot.org/uniprot}entry/{http://uniprot.org/uniprot}feature"):
        beginAttr = feature.find(
            "{http://uniprot.org/uniprot}location/{http://uniprot.org/uniprot}begin")
        if beginAttr is not None:
            begin = beginAttr.attrib['position']
            end = feature.find(
                "{http://uniprot.org/uniprot}location/{http://uniprot.org/uniprot}end").attrib['position']
            description = feature.attrib['description']
            rangesWithLabels[description] = ((begin, end))
    return rangesWithLabels


print(getRange())
