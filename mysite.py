from Bio.Seq import Seq
from Bio import SeqIO

class Site:
    def __init__(self, desc, seq):
        self.desc = desc    # description line
        self.seq = seq      # sequence

    def get_pos(self):
        """Returns the position of the site. NOT EXPANDED SEQUENCE"""
        try: return self.pos
        except AttributeError:
            wds = self.desc.split()
            self.pos = int(wds[wds.index('Abs_pos:')+1]) - 1    # start from 0
            return self.pos
    
    def get_strand(self):
        try: return self.strand
        except AttributeError:
            wds = self.desc.split()
            self.strand = int(wds[wds.index('Strand:')+1])
            return self.strand

    def get_orthologous_site(self):
        try: return self.orthologous_site
        except AttributeError:
            wds = self.desc.split()
            self.orthologous_site = wds[wds.index('RBH:')+1]
            return self.orthologous_site

    def get_region(self):
        try: return self.region
        except AttributeError:
            wds = self.desc.split()
            self.region = wds[wds.index('Region:')+1]
            return self.region

    def get_site_type(self):
        try: return self.site_type
        except AttributeError:
            wds = self.desc.split()
            self.site_type = wds[wds.index('Site_type:')+1]
            return self.site_type


def read_site_file(filename):
    """Read file of list of sites in fasta format. For each site, the file
    contains one line of description (which is parsed to get attributes) and
    sequence information.
    Return list of Site objects.
    """
    sites = []
    with open(filename) as f:
        while True:
            desc = f.readline().strip()
            if not desc: break # EOF
            assert desc.startswith(">")
            seq = f.readline().strip()
            new_site = Site(desc, seq)
            sites.append(new_site)
    return sites


