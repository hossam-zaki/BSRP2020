import urllib.parse
import urllib.request

url = 'https://www.uniprot.org/uniprot/'

params = {
    'query': 'gene:rad51b organism:human',
    'format': 'tab',
    'columns': 'id',
    'limit': 1
}
data = urllib.parse.urlencode(params)
print(data)
data = data.encode('utf-8')
req = urllib.request.Request(url, data)
with urllib.request.urlopen(req) as f:
    response = f.read()
print(response)
