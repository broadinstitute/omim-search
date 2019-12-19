"""
This script retrieves OMIM data from omim.org and parses/converts relevant fields into a tsv table.

==================
OMIM DATA SOURCES:
==================
OMIM provides data through an API (https://omim.org/help/api) and in downloadable files (https://omim.org/downloads/)
The geneMap API endpoint provides only gene symbols and not the Ensembl gene id, while
genemap2.txt provides both, so we use genemap2.txt as the data source.


API endpoints:
-------------
http://api.omim.org/api/geneMap?chromosome=1
   returns a list of 'geneMap' objects - each representing a
   mimNumber, geneSymbols, geneName, comments, geneInheritance, and a phenotypeMapList
   which contains one or more mimNumber, phenotypeMimNumber, phenotype description, and
   phenotypeInheritance

http://api.omim.org/api/entry?mimNumber=612367&format=json&include=all
   returns detailed info on a particular mim id

Files:
-----
mim2gene.txt - contains basic info on mim numbers and their relationships.

For example:
     100500  moved/removed
     100600  phenotype
     100640  gene    216     ALDH1A1 ENSG00000165092,ENST00000297785
     100650  gene/phenotype  217     ALDH2   ENSG00000111275,ENST00000261733

genemap2.txt - contains chrom, gene_start, gene_end, cyto_location, mim_number,
    gene_symbols, gene_name, approved_symbol, entrez_gene_id, ensembl_gene_id, comments, phenotypes,
    mouse_gene_id  -  where phenotypes contains 1 or more phenotypes in the form
    { description }, phenotype_mim_number (phenotype_mapping_key), inheritance_mode;

Example genemap2.txt record:

   # Chromosome    Genomic Position Start    Genomic Position End    Cyto Location    Computed Cyto Location    Mim Number    Gene Symbols    Gene Name    Approved Symbol    Entrez Gene ID    Ensembl Gene ID    Comments    Phenotypes    Mouse Gene Symbol/ID
   chr1    2019328    2030752    1p36.33        137163    GABRD, GEFSP5, EIG10, EJM7    Gamma-aminobutyric acid (GABA) A receptor, delta    GABRD    2563    ENSG00000187730        {Epilepsy, generalized, with febrile seizures plus, type 5, susceptibility to}, 613060 (3), Autosomal dominant; {Epilepsy, idiopathic generalized, 10}, 613060 (3), Autosomal dominant; {Epilepsy, juvenile myoclonic, susceptibility to}, 613060 (3), Autosomal dominant    Gabrd (MGI:95622)

"""
import argparse
import json
import logging
import os
import re
import requests
import datetime
from tqdm import tqdm
import urllib.request

logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s')
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

OMIM_ENTRIES_URL = 'https://api.omim.org/api/entry?apiKey={omim_key}&include=geneMap&format=json&mimNumber={mim_numbers}'

OMIM_PHENOTYPE_MAP_METHOD_CHOICES = {
    1: 'the disorder is placed on the map based on its association with a gene, but the underlying defect is not known.',
    2: 'the disorder has been placed on the map by linkage; no mutation has been found.',
    3: 'the molecular basis for the disorder is known; a mutation has been found in the gene.',
    4: 'a contiguous gene deletion or duplication syndrome, multiple genes are deleted or duplicated causing the phenotype.',
}


def get_remote_file_size(url):
    if url.startswith("http"):
        response = requests.head(url)
        return int(response.headers.get('Content-Length', '0'))
    elif url.startswith("ftp"):
        return 0  # file size not yet implemented for FTP
    else:
        raise ValueError("Invalid url: {}".format(url))


def download_file(url, to_dir=".", verbose=True):
    """Download the given file and returns its local path.
     Args:
        url (string): HTTP or FTP url
     Returns:
        string: local file path
    """

    if not (url and url.startswith(("http://", "https://", "ftp://"))):
        raise ValueError("Invalid url: {}".format(url))
    local_file_path = os.path.join(to_dir, os.path.basename(url))
    remote_file_size = get_remote_file_size(url)
    if os.path.isfile(local_file_path) and os.path.getsize(local_file_path) == remote_file_size:
        logger.info("Re-using {} previously downloaded from {}".format(local_file_path, url))
        return local_file_path

    logger.info(f"Downloading {url}. File size: {remote_file_size}")
    input_iter = urllib.request.urlopen(url)
    if verbose:
        logger.info("Downloading {} to {}".format(url, local_file_path))
        input_iter = tqdm(input_iter, unit=" data" if url.endswith("gz") else " lines")

    with open(local_file_path, 'w') as f:
        f.writelines(input_iter)

    input_iter.close()

    return local_file_path


CHROMOSOMES = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y', 'M']
CHROM_TO_CHROM_NUMBER = {chrom: i for i, chrom in enumerate(CHROMOSOMES)}


def get_xpos(chrom, pos):
    """Compute single number representing this chromosome and position.

    Args:
        chrom (string): examples '1', 'Y', 'M'
        pos (integer): genomic position on chromosome
    """
    if chrom not in CHROM_TO_CHROM_NUMBER:
        chrom = chrom.replace('chr', '')
        if chrom.startswith('M'):
            chrom = 'M'
        if chrom not in CHROM_TO_CHROM_NUMBER:
            raise ValueError("Invalid chromosome: %s" % (chrom,))

    if pos < 1 or pos > 3e8:
        raise ValueError("Invalid position: %s" % (pos,))

    return (1 + CHROM_TO_CHROM_NUMBER[chrom])*int(1e9) + pos


def get_file_header(f):
    header_fields = None
    for i, line in enumerate(f):
        line = line.rstrip('\r\n')
        if line.startswith("# Chrom") and header_fields is None:
            header_fields = [c.lower().replace(' ', '_') for c in line.lstrip('# ').split('\t')]
            break
        elif not line or line.startswith("#"):
            continue
        elif line.startswith('This account is inactive') or line.startswith('This account has expired'):
            raise Exception(line)
        elif header_fields is None:
            raise ValueError("Header row not found in genemap2 file before line {}: {}".format(i, line))

    return header_fields


def parse_genemap2_records(omim_line_fields):
    # skip commented rows
    if len(omim_line_fields) == 1:
        yield None

    else:
        # rename some of the fields
        output_record = {}
        output_record['chrom'] = omim_line_fields['chromosome'].replace("chr", "")
        output_record['start'] = int(omim_line_fields['genomic_position_start']) or 1  # 'or 1' replaces pos=0 with 1
        output_record['end'] = int(omim_line_fields['genomic_position_end']) or 1
        output_record['cyto'] = omim_line_fields['cyto_location']
        output_record['gene_id'] = omim_line_fields['ensembl_gene_id']
        output_record['mim_number'] = int(omim_line_fields['mim_number'])
        output_record['gene_symbols'] = ", ".join([s for s in [omim_line_fields['approved_symbol'].strip()] + list(omim_line_fields['gene_symbols'].split(",")) if s])
        output_record['gene_description'] = omim_line_fields['gene_name']
        output_record['comments'] = omim_line_fields['comments']
        output_record['mouse_gene_id'] = omim_line_fields['mouse_gene_symbol/id']

        phenotype_field = omim_line_fields['phenotypes'].strip()

        record_with_phenotype = None
        for phenotype_match in re.finditer("[\[{ ]*(.+?)[ }\]]*(, (\d{4,}))? \(([1-4])\)(, ([^;]+))?;?", phenotype_field):
            # Phenotypes example: "Langer mesomelic dysplasia, 249700 (3), Autosomal recessive; Leri-Weill dyschondrosteosis, 127300 (3), Autosomal dominant"

            record_with_phenotype = dict(output_record)  # copy
            record_with_phenotype["phenotype_description"] = phenotype_match.group(1)
            record_with_phenotype["phenotype_mim_number"] = int(phenotype_match.group(3)) if phenotype_match.group(
                3) else None
            record_with_phenotype["phenotype_map_method"] = phenotype_match.group(4)
            record_with_phenotype["phenotype_inheritance"] = phenotype_match.group(6) or None

            # basic checks
            if len(record_with_phenotype["phenotype_description"].strip()) == 0:
                raise ValueError("Empty phenotype description: {}".format(json.dumps(omim_line_fields)))

            if int(record_with_phenotype["phenotype_map_method"]) not in OMIM_PHENOTYPE_MAP_METHOD_CHOICES:
                raise ValueError("Unexpected value (%s) for phenotype_map_method: %s" % (
                    record_with_phenotype["phenotype_map_method"], phenotype_field))

            yield record_with_phenotype

        if record_with_phenotype is None:
            if len(phenotype_field) > 0:
                raise ValueError("No phenotypes found: {}".format(json.dumps(omim_line_fields)))
            else:
                yield output_record


def add_phenotypic_series(records, omim_key):
    logger.info('Adding phenotypic series information via the OMIM API.')
    mim_numbers = set()
    for omim_record in records:
        if omim_record.get("phenotype_mim_number"):
            mim_numbers.add(str(omim_record["mim_number"]))
    mim_numbers = list(mim_numbers)

    mim_number_to_phenotypic_series = {}
    for i in tqdm(range(0, len(mim_numbers), 20), unit=" batches"):
        logger.debug('Fetching entries {}-{}'.format(i, i + 20))
        entries_to_fetch = mim_numbers[i:i + 20]
        response = requests.get(OMIM_ENTRIES_URL.format(omim_key=omim_key, mim_numbers=','.join(entries_to_fetch)))
        if not response.ok:
            raise Exception('Request failed with {}: {}'.format(response.status_code, response.reason))

        entries = response.json()['omim']['entryList']
        if len(entries) != len(entries_to_fetch):
            raise Exception('Expected {} omim entries but recieved {}'.format(len(entries_to_fetch), len(entries)))

        for entry in entries:
            mim_number = entry['entry']['mimNumber']
            for phenotype in entry['entry'].get('geneMap', {}).get('phenotypeMapList', []):
                phenotypic_series_number = phenotype['phenotypeMap'].get('phenotypicSeriesNumber')
                if phenotypic_series_number:
                    mim_number_to_phenotypic_series[mim_number] = phenotypic_series_number

    for omim_record in records:
        omim_record["phenotypic_series_number"] = mim_number_to_phenotypic_series.get(omim_record["mim_number"])

    logger.info(f'Found {len(mim_number_to_phenotypic_series)} records with phenotypic series')


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--omim-key', help="OMIM key provided with registration", default=os.environ.get("OMIM_KEY"))
    args = p.parse_args()

    if not args.omim_key:
        p.error("--omim-key not specified and OMIM_KEY env var not set.")
    return args


def main():
    args = parse_args()

    try:
        file_path = download_file(f"https://data.omim.org/downloads/{args.omim_key}/genemap2.txt")
    except Exception as e:
        logger.error(e)
        file_path = "genemap2.txt"


    file_mod_timestamp = datetime.datetime.fromtimestamp(os.path.getmtime(file_path))

    logger.info(f'Parsing {file_path}')
    records = []
    with open(file_path) as f:
        header_fields = get_file_header(f)

        for line in tqdm(f, unit=" records"):
            omim_line_fields = dict(zip(header_fields, line.rstrip('\r\n').split('\t')))
            for record in parse_genemap2_records(omim_line_fields):
                if record is None:
                    continue

                records.append(record)

    add_phenotypic_series(records, args.omim_key)

    data_matrix = []
    for record in records:
        record['xstart'] = get_xpos(record['chrom'], record['start'])
        record['xend'] = get_xpos(record['chrom'], record['end'])
        record['locus'] = f"{record['chrom']}:{record['start']}-{record['end']}"

        matrix_row = []
        for key in [
            'xstart',
            'xend',
            'locus',
            'cyto',
            'gene_symbols',
            'gene_id',
            'mim_number',
            'phenotype_mim_number',
            'phenotype_inheritance',
            'phenotype_map_method',
            'phenotypic_series_number',
            'gene_description',
            'phenotype_description',
            'comments',
            'mouse_gene_id',
        ]:
            value = record.get(key)
            matrix_row.append(value if value is not None else '')

        data_matrix.append(matrix_row)

    json_string = json.dumps({
        "data": data_matrix,
        "createdDate": file_mod_timestamp.strftime('%Y-%m-%d %H:%M:%S %p')
    })

    with open("omim.json", "wt") as f:
        f.write(json_string)

    logger.info(f"Generated omim.json with {len(records)} records")


if __name__ == "__main__":
    main()
