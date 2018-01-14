import pandas as pd
import os


class UseCosmic:
    def __init__(self):
        self.cosmicdf = self.load_df()

    def get_cosmic(self):
        """get cosmic df locatiaon"""
        user_input = input("No cosmic file designed. Enter the path of cosmic your file: ")
        cosmic_file = user_input
        cosmic_file = str(cosmic_file)
        path = os.path.dirname(__file__)
        file = open('{}/config.py'.format(path), 'w')
        file.write('cosmic_file = "{}"'.format(cosmic_file))
        file.close()
        return cosmic_file

    def compare(self, genelist):
        """Compare genes in a list to df list"""
        cosmiclist = self.cosmicdf['Gene Symbol'].tolist()
        return frozenset(genelist).intersection(cosmiclist)

    def custom_compare(self, genelist, genes_of_intrest):
        """ Compare a random gene list NOT CURRENTLY USED"""
        return frozenset(genelist).intersection(genes_of_intrest)

    def load_df(self):
        """Loads cosmic_file from config.py or runs get_cosmic else  """
        try:
            from .config import cosmic_file
            cosmicdf = pd.read_csv(cosmic_file)
        except ImportError:
            cosmic_file = self.get_cosmic()
            cosmicdf = pd.read_csv('{}'.format(cosmic_file))
        return cosmicdf



