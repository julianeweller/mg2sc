import pysam
import regex as re
import collections.abc
from scipy.sparse import csr_matrix
from scipy.io import mmwrite
import csv
import logging
import os
    

def mg2sc(bamfile, mgfile, dbfile, outdir):
    """ Main Function. 
    Creates a sparse matrix with transcript count per organism for each cell."""

    # Generate variables based on input
    matrixfile = outdir + 'matrix.mtx'
    cellfile = outdir + 'barcodes.tsv'
    taxfile = outdir + 'taxids.tsv'
    dbfile = os.path.join(dbfile, 'inspect.txt')
    dbfile_out = outdir + 'inspect_db.txt'

    # Extract taxonomy IDs for each transcript
    mg_dict = extract_ids(bamfile, mgfile)

    # Find most frequent taxonomy for each transcript
    map_nested_dicts(mg_dict, most_frequent)

    # Make sparse matrix
    rows, cols, vals, cell_list, taxid_list = dict2lists(twist_dict(mg_dict))
    sparsematrix = csr_matrix((vals, (rows, cols)))

    # Get ncbi name for taxonomy ID
    taxdict = krakenID2dict(dbfile, taxid_list)
    taxname_list = [taxdict[k] for k in taxid_list]

    # store sparse matrix
    mmwrite(matrixfile, sparsematrix)
    
    # Store list of cell barcodes
    with open(cellfile, 'w') as f_output:
        tsv_output = csv.writer(f_output, delimiter='\n')
        tsv_output.writerow(cell_list)
    
    # Store list of taxonomy IDs
    data = zip(taxid_list, taxname_list)
    with open(taxfile, 'w') as f_output:
        tsv_output = csv.writer(f_output, delimiter='\t')
        for idx, tax  in data:
            tsv_output.writerow([idx, tax])
    
    # Store reference database hierarchy
    with open(dbfile) as f:
        with open(dbfile_out, "w") as f1:
            for line in f:
                f1.write(line) 

def extract_ids(bamfile, krakenfile): 
    """
    Builds a nested dictionary with KRAKEN2 taxonomy code for each transcript and the cell it belongs to.
    Input:  Output from KRAKEN2, .bam file with unmapped reads
    Output: {cellbarcode: {transcriptbarcode: krakentaxonomyID}}
    """
    line = 0
    skipped = 0
    # Store extracted information in nested dictionary {cellbarcode:{transcriptbarcode: taxonomyID}}
    nested_dict = {}
    
    # Iterate simultanously through bam and kraken file
    for sread,kread in zip(pysam.AlignmentFile(bamfile, "rb"),open(krakenfile,"r")):
        
        # count the total number of reads analysed
        line += 1
        
        # Check that read names in kraken and bam file match
        if sread.query_name != kread.split('\t')[1]:
            skipped += 1
            logging.warning("sam file read name and metagenomicsfile read name don't match and are therefore excluded: sam: {}, kraken: {}".format(sread.query_name, kread.split('\t')[1]))
            continue

        # Get cell barcode and UMI from bam file
        try:
            sread_CB = sread.get_tag('CB')
            sread_UB = sread.get_tag('UB')
        except:
            # some reads don't have a cellbarcode or transcript barcode. They can be skipped.
            skipped += 1
            continue
            
        # Get taxonomy ID from kraken file
        kread_taxid = kread.split('\t')[2]
        if (type(kread_taxid) != int) and (kread_taxid.isdigit() == False):
            try:
                # sometimes, the taxonomy is name (taxid #), sometimes it's just the number
                kread_taxid = re.search('\(([^)]+)', kread_taxid).group(1)[6:]
            except:
                # in this case, something is wrong!
                logging.debug("Here is an error. TaxID: {}".format(kread_taxid))
                sys.exit()

        # Make nested dictionary with cells and transcripts
        if sread_CB in nested_dict:
            # If cell and transcript exist, add taxonomy ID to list
            if sread_UB in nested_dict[sread_CB]:
                nested_dict[sread_CB][sread_UB].append(kread_taxid)
            # Otherwise create transcript dictionary for cell
            else:
                nested_dict[sread_CB][sread_UB] = [kread_taxid]
        else:
            # if cell doesn't exist, create cell and transcript dictionary with kraken id
            nested_dict[sread_CB] = {sread_UB: [kread_taxid]}
    
    # Output control values
    logging.info("total reads: {}, skipped reads: {}".format(line,skipped))
    
    return nested_dict

def most_frequent(List):
    """Finds the most frequent element in a list"""
    return max(set(List), key = List.count)

def map_nested_dicts(ob, func):
    """ Applys a map to the inner item of nested dictionaries """
    for k, v in ob.items():
        if isinstance(v, collections.abc.Mapping):
            map_nested_dicts(v, func)
        else:
            ob[k] = func(v)

def twist_dict(nested):
    """ Make count dictionary with {cellbarcode : {taxonomyID : transcriptcount}} """
    newdict = {}
    for ckey, tdict in nested.items():
        for tkey, kvalue in tdict.items():
            if ckey in newdict:
                if kvalue in newdict[ckey]:
                    newdict[ckey][kvalue] += 1
                else:
                    newdict[ckey][kvalue] = 1
            else:
                newdict[ckey] = {kvalue: 1}
    return(newdict)

def dict2lists(nested):
    """ Returns lists for sparse matrix """
    rows = [] # cell coordinate
    columns = [] # taxonomy id coordinate
    values = [] # count

    cell_list = [] # same order as rows
    taxid_list = [] # same order as columns

    j = 0

    for ckey, taxdict in nested.items():
        for taxkey, count in taxdict.items():
            try:
                k = taxid_list.index(taxkey)
            except:
                taxid_list.append(taxkey)
                k = taxid_list.index(taxkey)
                
            rows.append(k)
            columns.append(j)
            values.append(count) 
            
        # increase cell coordinate by 1
        cell_list.append(ckey)
        j += 1
    
    return rows, columns, values, cell_list, taxid_list

def krakenID2dict(dbfile, taxid_list):
    """
    Get name for each taxonomy ID from kraken database
    """
    # iterate through inspect file and lookup taxonomy ids
    k=0
    taxdict = {'0': 'unclassified'}

    with open(dbfile) as f:
        for line in f:
            if line.startswith("#"):
                continue
            
            # string to list
            line = line[:-1].split('\t')
            taxid_db = line[4]
            taxname = line[5].lstrip()
            
            if taxid_db in taxid_list:
                taxdict[taxid_db] = taxname
    
    return taxdict

def extract_taxref(file):
    """ 
    Extract taxonomy reference for each read.
    Input:  viral track output .bam file
    Output: dictionary with {readname: taxonomy ID}, list of unique taxonomy IDs
    """
    # extract taxref for each read
    tdict = {}
    line = 0
    skipped = 0
    taxref_list = set('0')
    
    for read in pysam.AlignmentFile(file, "rb"):
        # count the total number of reads analysed
        line += 1
        try:
            # Extract readname and taxonomy reference
            taxref = read.to_dict().get('ref_name').split('|')[1]
            taxref_list.add(taxref)
            tdict[read.query_name] = taxref
        except:
            # in case some reads are unmapped or don't work
            skipped += 1
    logging.info("Reads in ViralTrack output: {}, reads without taxonomy reference or that failed: {}".format(line, skipped))
    return(tdict, taxref_list)

def extract_bc(file):
    """ 
    Extracts cellbarcode and UMI for each readname
    Input:  unmapped .bam file
    Output: dictionary with {readname: [cellbarcode, UMI]}
    """
    # extract UB and CB for each read
    bcdict = {}
    line = 0
    skipped = 0

    for read in pysam.AlignmentFile(file, "rb"):
        # count the total number of reads analysed
        line += 1
        # Get cell barcode and UMI from bam file
        try:
            # Extract readname, cell barcode and UMI
            bcdict[read.query_name] = [read.get_tag('CB'),read.get_tag('UB')]
        except:
            # some reads don't have a cellbarcode or transcript barcode. They can be skipped.
            skipped += 1
            continue

    logging.info("Reads in original bam file: {}, reads without cellbarcode or UMI: {}".format(line, skipped))
    return(bcdict)

# if __name__ == '__main__':
#     # add here argparse for files
#     mgfile =''
#     bamfile = ''
#     dbfile = ''
#     outdir = ''
#     mg2sc(bamfile, mgfile, dbfile, outdir)