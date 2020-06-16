# Copyright (c) 2017 The Ontario Institute for Cancer Research. All rights
# reserved.
#
# This program and the accompanying materials are made available under the
# terms of the GNU Public License v3.0.
# You should have received a copy of the GNU General Public License along
# with this program. If not, see <http://www.gnu.org/licenses/>.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
# IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
# THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED.
# IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY
# DIRECT,  INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
# HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
# IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
# OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
# ADVISED OF THE  POSSIBILITY OF SUCH DAMAGE.
"""
query.py
This script demonstrates running a simple PQL query against the ICGC data
portal with the icgc module.
"""
from __future__ import absolute_import, print_function

import json
import urllib

import icgc
import matplotlib.pyplot as plt
import pandas as pd


def run():
    """
    Demonstrate PQL by displaying 1 of each request type as JSON output
    """
    numberOfSeqs = json.load(urllib.request.urlopen(
        "http://dcc.icgc.org/api/v1/genes/ENSG00000182185/mutations/count"))
    print(numberOfSeqs)
    counter = 0
    to_nearest_hunderd = 101 - (numberOfSeqs % 100)
    placeInGenome = []
    numberOfOccurences = []
    for i in range(0, numberOfSeqs+to_nearest_hunderd, 100):
        try:
            json_file = json.load(urllib.request.urlopen(
                f"https://dcc.icgc.org/api/v1/genes/ENSG00000182185/mutations?from={i}&size=100"))
        except:
            print(f"{i} got messed")
        print(i)
        for hit in json_file['hits']:
            if hit['type'] == 'single base substitution':
                placeInGenome.append(hit['start'])
                numberOfOccurences.append(hit['affectedDonorCountTotal'])
    plt.stem(placeInGenome, numberOfOccurences)
    plt.title("RAD51B Mutation Map")
    plt.xlabel('Position')
    plt.ylabel('# of Donors')

    plt.show()

    # response = icgc.query(request_type='genes',
    #                       pql='')
    # print(response)
    # json_file = json.load(urllib.request.urlopen(
    #     "http://dcc.icgc.org/api/v1/genes/ENSG00000182185/mutations?filters=%7B%7D&from=1&size=10&sort=affectedDonorCountFiltered&order=desc"))
    # for doc in json_file['hits']:
    #     print(doc['id'])
    #     # if ">" in doc['mutation']:
    #     #     print("yee")


if __name__ == '__main__':
    run()
