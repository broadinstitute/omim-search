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

http://api.omim.org/api/entry?apiKey={omim_key}&include=geneMap&format=json&mimNumber={mim_numbers}
http://api.omim.org/api/entry?apiKey={omim_key}&include=all&format=json&mimNumber=612367
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
import collections
import datetime
import json
import logging
import os
import re
import requests

from tqdm import tqdm
import urllib.request

from pyliftover.liftover import LiftOver

#import hail as hl
#hl.get_reference('GRCh38').add_liftover("gs://hail-common/references/grch38_to_grch37.over.chain.gz", 'GRCh37')


logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s')
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

OMIM_ENTRIES_URL = 'https://api.omim.org/api/entry?apiKey={omim_key}&include=all&format=json&mimNumber={mim_numbers}'

OMIM_PHENOTYPE_MAP_METHOD_CHOICES = {
    1: 'the disorder is placed on the map based on its association with a gene, but the underlying defect is not known.',
    2: 'the disorder has been placed on the map by linkage; no mutation has been found.',
    3: 'the molecular basis for the disorder is known; a mutation has been found in the gene.',
    4: 'a contiguous gene deletion or duplication syndrome, multiple genes are deleted or duplicated causing the phenotype.',
}

LIFTOVER_HG38_TO_HG19 = LiftOver('hg38', 'hg19')


def download_file(url, to_dir=".", force_download=False, verbose=True):
    """Download the given file and returns its local path.
     Args:
        url (string): HTTP or FTP url
     Returns:
        string: local file path
    """

    if not (url and url.startswith(("http://", "https://", "ftp://"))):
        raise ValueError("Invalid url: {}".format(url))

    local_file_path = os.path.join(to_dir, os.path.basename(url))
    if not force_download and os.path.isfile(local_file_path) and  os.path.getsize(local_file_path) > 100:
        return local_file_path

    logger.info(f"Downloading {url} to {local_file_path}")

    input_iter = urllib.request.urlopen(url)
    if verbose:
        input_iter = tqdm(input_iter, unit=" data" if url.endswith("gz") else " lines")

    with open(local_file_path, 'w') as f:
        f.writelines((line.decode('utf-8') for line in input_iter))

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
        output_record['gene_symbols'] = ", ".join(sorted(set([s for s in [omim_line_fields['approved_symbol'].strip()] + list(omim_line_fields['gene_symbols'].split(",")) if s])))
        output_record['gene_description'] = omim_line_fields['gene_name']
        output_record['comments'] = omim_line_fields['comments']
        output_record['mouse_gene_id'] = omim_line_fields['mouse_gene_symbol/id']

        phenotype_field = omim_line_fields['phenotypes'].strip()

        record_with_phenotype = None
        for phenotype_match in re.finditer("[\[{ ]*(.+?)[ }\]]*(, (\d{4,}))? \(([1-4])\)(, ([^;]+))?;?", phenotype_field):
            # Phenotypes example: "Langer mesomelic dysplasia, 249700 (3), Autosomal recessive; Leri-Weill dyschondrosteosis, 127300 (3), Autosomal dominant"

            record_with_phenotype = dict(output_record)  # copy
            record_with_phenotype["phenotype_description"] = phenotype_match.group(1)
            record_with_phenotype["phenotype_mim_number"] = int(phenotype_match.group(3)) if phenotype_match.group(3) else None
            record_with_phenotype["phenotype_map_method"] = int(phenotype_match.group(4))
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
    logger.info('Adding phenotypic_series_number, text, date_created and date_updated columns via the OMIM API.')
    mim_numbers = set()
    for omim_record in records:
        if omim_record.get("mim_number"):
            mim_numbers.add(str(omim_record["mim_number"]))
    mim_numbers = list(mim_numbers)

    # query the API for entries for each MIM number
    phenotypic_series_found_counter = 0
    mim_number_to_api_info = collections.defaultdict(dict)
    for i in tqdm(range(0, len(mim_numbers), 20), unit=" batches"):
        logger.debug('Fetching entries {}-{}'.format(i, i + 20))
        entries_to_fetch = mim_numbers[i:i + 20]
        response = requests.get(OMIM_ENTRIES_URL.format(omim_key=omim_key, mim_numbers=','.join(entries_to_fetch)))
        if not response.ok:
            raise Exception('Request failed with {}: {}'.format(response.status_code, response.reason))

        entries = response.json()['omim']['entryList']
        if len(entries) != len(entries_to_fetch):
            raise Exception('Expected {} omim entries but received {}'.format(len(entries_to_fetch), len(entries)))

        for entry in entries:
            mim_number = entry['entry']['mimNumber']

            datetime_created = datetime.datetime.fromtimestamp(int(entry['entry']['epochCreated']))
            mim_number_to_api_info[mim_number]['date_created'] = datetime_created.strftime("%Y-%m-%d")

            datetime_updated = datetime.datetime.fromtimestamp(int(entry['entry']['epochUpdated']))
            mim_number_to_api_info[mim_number]['date_updated'] = datetime_updated.strftime("%Y-%m-%d")

            for phenotype in entry['entry'].get('geneMap', {}).get('phenotypeMapList', []):
                phenotypic_series_number = phenotype['phenotypeMap'].get('phenotypicSeriesNumber')
                if phenotypic_series_number:
                    mim_number_to_api_info[mim_number]['phenotypic_series_number'] = phenotypic_series_number
                    phenotypic_series_found_counter += 1
                    break

            mim_number_to_api_info[mim_number]['text'] = ""
            for textSection in entry['entry'].get('textSectionList', []):
                if textSection.get('textSection', {}).get('textSectionName') in ("text", "description"):
                    mim_number_to_api_info[mim_number]['text'] += textSection.get('textSection', {}).get('textSectionContent')

    # transfer API info to the omim_record for each mim number
    for omim_record in records:
        omim_record.update(mim_number_to_api_info.get(omim_record["mim_number"], {}))

    logger.info(f'Found {phenotypic_series_found_counter} records with phenotypic series')


def convert_hg38_to_h19_using_pyliftover(hg38_chrom, hg38_start, hg38_end):
    hg19_start_coord = LIFTOVER_HG38_TO_HG19.convert_coordinate(hg38_chrom, int(hg38_start))
    if not hg19_start_coord or not hg19_start_coord[0]:
        return None, None, None

    hg19_end_coord = LIFTOVER_HG38_TO_HG19.convert_coordinate(hg38_chrom, int(hg38_end))
    if not hg19_end_coord or not hg19_end_coord[0]:
        return None, None, None

    hg19_chrom = hg19_start_coord[0][0].lstrip('chr')
    if hg19_chrom not in CHROMOSOMES:
        return None, None, None

    hg19_start = min(hg19_start_coord[0][1], hg19_end_coord[0][1])
    hg19_end = max(hg19_start_coord[0][1], hg19_end_coord[0][1])

    return hg19_chrom, hg19_start, hg19_end


def get_hg19_coordinates_for_gene_id(gene_id): # A simple function to use requests.post to make the API call. Note the json= section.
    request = requests.post('https://gnomad.broadinstitute.org/api', json={
        'query': f'{{gene(gene_id:"{gene_id}",reference_genome:GRCh37){{chrom, start, stop}}}}'
    })
    if request.status_code == 200:
        results = request.json()
        gene = results.get('data', {}).get('gene', {}) or {}
        return gene.get('chrom'), gene.get('start'), gene.get('stop')


def convert_hg38_to_h19(record):
    hg38_chrom = 'chr{}'.format(record['chrom'].lstrip('chr'))
    hg38_start = int(record['start'])
    hg38_end = int(record['end'])

    hg19_chrom, hg19_start, hg19_end = convert_hg38_to_h19_using_pyliftover(hg38_chrom, hg38_start, hg38_end)
    if hg19_chrom is None and record['gene_id']:
        hg19_chrom, hg19_start, hg19_end = get_hg19_coordinates_for_gene_id(record['gene_id'])

    return hg19_chrom, hg19_start, hg19_end


def add_hg19_coords(records):
    failed_liftover = []
    failed_liftover_gene_ids = []
    for record in tqdm(records, unit=" records"):

        hg19_chrom, hg19_start, hg19_end = convert_hg38_to_h19(record)

        record["liftover_to_hg19_failed"] = False
        if hg19_chrom is None:
            hg19_chrom, hg19_start, hg19_end = record['chrom'], int(record['start']), int(record['end'])

            #logger.info(f"Couldn't lift {record['locus']} from hg38 => hg19. Gene ID: {record['gene_id']}")
            failed_liftover.append(record)
            record["liftover_to_hg19_failed"] = True
            if record['gene_id']:
                failed_liftover_gene_ids.append(record['gene_id'])

        record['xstart_hg19'] = get_xpos(hg19_chrom, hg19_start)
        record['xend_hg19'] = get_xpos(hg19_chrom, hg19_end)
        record['locus_hg19'] = f"{hg19_chrom}:{hg19_start}-{hg19_end}"
        if record['liftover_to_hg19_failed']:
            record['locus_hg19'] += " ** failed liftover from hg38 to hg19. Reusing hg38 coords."

    logger.info(f"Liftover failed for {len(failed_liftover)} out of {len(records)} records ({100*len(failed_liftover)/len(records):0.1f}%)")
    logger.info(f"{len(set(failed_liftover_gene_ids))} had gene ids: {set(failed_liftover_gene_ids)}")


"""
def add_hg19_coords_using_hail(records):
    for record in tqdm(records, unit=" records"):
        hg19_coords = hl.eval(hl.liftover(
            hl.locus_interval(f"chr{record['chrom']}", record['start'], record['end'], True, True, 'GRCh38'),
            'GRCh37',
            include_strand=True,
        ))

        if hg19_coords is not None:
            #strand = "-" if hg19_coords.is_negative_strand else "+"
            record['locus_hg19'] = f"{hg19_coords.result.start.contig}:{hg19_coords.result.start.position}-{hg19_coords.result.end.position}"
        else:
            record['locus_hg19'] = None
"""


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--omim-key', help="OMIM key provided with registration", default=os.environ.get("OMIM_KEY"))
    p.add_argument('--force', help="Download from OMIM even if local genemap2 file already exists", action="store_true")
    args = p.parse_args()

    if not args.omim_key:
        p.error("--omim-key not specified and OMIM_KEY env var not set.")
    return args


def main():
    args = parse_args()

    try:
        file_path = download_file(f"https://data.omim.org/downloads/{args.omim_key}/genemap2.txt", force_download=args.force)
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

    for record in records:
        record['xstart'] = get_xpos(record['chrom'], record['start'])
        record['xend'] = get_xpos(record['chrom'], record['end'])
        record['locus'] = f"{record['chrom']}:{record['start']}-{record['end']}"

    add_hg19_coords(records)
    add_phenotypic_series(records, args.omim_key)

    data_matrix = []
    for record in records:
        matrix_row = []
        for key in [
            'mim_number',                   # 0
            'phenotype_mim_number',         # 1
            'phenotypic_series_number',     # 2
            'phenotype_inheritance',        # 3

            'locus',           # 4
            'locus_hg19',      # 5

            'cyto',            # 6
            'gene_symbols',    # 7
            'gene_id',         # 8
            'gene_description',        # 9
            'phenotype_description',   # 10
            'date_created',            # 11
            'date_updated',            # 12
            'mouse_gene_id',           # 13
            'text',                    # 14
            'comments',                # 15

            'xstart',
            'xend',
            'xstart_hg19',
            'xend_hg19',
            'phenotype_map_method',
            #'liftover_to_hg19_failed',
        ]:
            value = record.get(key)
            matrix_row.append(value if value is not None else '')

        data_matrix.append(matrix_row)

    json_string = json.dumps({
        "data": data_matrix,
        "date_downloaded_from_omim": file_mod_timestamp.strftime('%Y-%m-%d %I:%M:%S %p')
    })

    with open("omim.json", "wt") as f:
        f.write(json_string)

    logger.info(f"Generated omim.json with {len(records)} records")


if __name__ == "__main__":
    main()
