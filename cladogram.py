from collections import deque
import paleobiodb_interface as pbdb
from paleobiodb_interface import clad_vocab as cv

def cladogram_lengths(taxon_name):
    # Initial query to get cladogram data
    taxa = pbdb.query_cladogram(taxon_name)

    # Load cladogram into stack LIFO (root on top)
    backlog = deque()
    backlog.extend(taxa)
    cladogram = {}
    while True:
        if len(backlog) <= 0:
            break
        taxon = backlog.popleft()
        if cv.FLAG in taxon and taxon[cv.FLAG] == 'B':
            cladogram[taxon[cv.ID]] = 0
        elif cv.PARENT_ID not in taxon:
            continue
        elif taxon[cv.PARENT_ID] not in cladogram.keys(): 
            backlog.append(taxon)
        else:
            cladogram[taxon[cv.ID]] = cladogram[taxon[cv.PARENT_ID]] + 1

    return cladogram

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Get cladogram lengths for a given taxon name from Paleobiology Database')
    parser.add_argument('taxon_name', type=str, help='Taxon name to query')
    args = parser.parse_args()

    lengths = cladogram_lengths(args.taxon_name)
    for taxon_id, length in lengths.items():
        print(f'Taxon ID: {taxon_id}, Length from root: {length}')