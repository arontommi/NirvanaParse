import argparse
import ijson
import pandas as pd
from Cosmic_comparison import cclookout


def parser():
    """Get arguments"""
    parserf = argparse.ArgumentParser()
    parserf.add_argument("-i", help="input file")
    parserf.add_argument("-o", help="output file")
    parserf.add_argument("-s", help="sample name")

    args = parserf.parse_args()
    return args


def get_coverage(dfcolumn, keep):
    """Get Coverage data """
    lol = []
    for line in dfcolumn:
        sv_ids = []
        try:
            for ii in line:
                try:
                    sv_ids.append(ii['{}'.format(keep)])
                except KeyError:
                    sv_ids.append([0, 0])
        except TypeError:
            sv_ids.append('[]')
        lol.append(sv_ids)
    return lol


def fix_coverage(dfcolumn):
    nref, nalt, tref, talt = [], [], [], []
    for row in dfcolumn:
        nref.append(list(row[0])[0])
        nalt.append(list(row[0])[1])
        tref.append(list(row[1])[0])
        talt.append(list(row[1])[1])
    return nref, nalt, tref, talt


def get_gene_stuff(column, ids):
    """Get Gene info"""
    uc = cclookout.UseCosmic()
    listoflists = []
    cancergenes = []
    for row in column:
        idlist = []
        try:
            for ii in row:
                try:
                    idlist.append(ii['{}'.format(ids)])
                except KeyError:
                    idlist.append('')
            idlist = [item for sublist in idlist for item in sublist]
            if len(idlist) <= 10:
                listoflists.append(idlist)
                cancergenes.append(uc.compare(idlist))
            else:
                cancergenes.append(uc.compare(idlist))
                listoflists.append(['To many values to keep'])
        except TypeError:
            listoflists.append([]) 
    return listoflists, cancergenes


def get_dict_stuff(column, ids):
    """Get all the other dict stuff"""
    listoflists = []
    for row in column:
        idlist = []
        try:
            for ii in row:
                try:
                    idlist.append(ii['{}'.format(ids)])
                except KeyError:
                    idlist.append('')
            if any(isinstance(sublist, list) for sublist in idlist):
                idlist = [item for sublist in idlist for item in sublist]
            if len(idlist) <= 10:
                listoflists.append(idlist)
            else:
                listoflists.append(['To many values to keep'])
        except TypeError:
            listoflists.append([]) 
    return listoflists


def getdata(filename):
    """Parse data to lists"""
    data = []
    print('Loading json output to a pandas dataframe')
    with open(filename, 'r') as f:
        objects = ijson.items(f, 'positions')
        for i in objects:
            data.append(i)
    print('Finish loading, starting filtering ')
    return data


def getmaindf(data):
    """create a df and transpose from data"""
    maindf = pd.DataFrame()
    for i in data[0]:
        df = pd.DataFrame.from_dict(i, orient="index")
        dft = df.T
        maindf = maindf.append(dft)
    return maindf


def main():
    """ run all the stuff """
    args = parser()
    filename = args.i
    data = getdata(filename)
    maindf = getmaindf(data)

    maindf['filters'] = maindf.filters.apply(', '.join)
    maindf = maindf[maindf['filters'] == 'PASS']
    maindf['sample'] = args.s

    maindf['Genes'], maindf['Cancer_Genes_in_region'] = get_gene_stuff(maindf['variants'], 'overlappingGenes')

    maindf['transcripts'] = get_dict_stuff(maindf['variants'], 'transcripts')
    maindf['Coverage N/T Paired'] = get_coverage(maindf['samples'], 'pairedEndReadCounts')
    maindf['Paired Normal ref depth'], maindf['Paired Normal alt depth'], maindf['Paired Tumor ref depth'],\
    maindf['Paired Tumor alt depth'] = fix_coverage(maindf['Coverage N/T Paired'])

    maindf['Coverage N/T Split'] = get_coverage(maindf['samples'], 'splitReadCounts')
    maindf['Split Normal ref depth'], maindf['Split Normal alt depth'], maindf['Split Tumor ref depth'],\
    maindf['Split Tumor alt depth'] = fix_coverage(maindf['Coverage N/T Split'])

    maindf['altAlleles'] = maindf.altAlleles.apply(', '.join)
    maindf['Genes'] = maindf.Genes.apply(', '.join)
    maindf['Cancer_Genes_in_region'] = maindf.Cancer_Genes_in_region.apply(', '.join)

    maindf['paired T frequency'] = (maindf['Paired Tumor alt depth'] / maindf['Paired Tumor ref depth']).round(2)
    maindf['split T frequency'] = (maindf['Split Tumor alt depth'] / maindf['Split Tumor ref depth']).round(2)

    keeps = ['sample',
             'chromosome',
             'position',
             'refAllele',
             'altAlleles',
             'cytogeneticBand',
             'Genes',
             'Cancer_Genes_in_region',
             'svLength',
             'jointSomaticNormalQuality',
             'paired T frequency',
             'split T frequency',
             'Paired Tumor ref depth',
             'Paired Tumor alt depth',
             'Split Tumor ref depth',
             'Split Tumor alt depth',
             'Paired Normal ref depth',
             'Paired Normal alt depth',
             'Split Normal ref depth',
             'Split Normal alt depth',
             'filters']

    maindf = maindf[keeps]
    maindf.to_csv('{}'.format(args.o), sep='\t', index=None)


if __name__ == '__main__':
    main()
