class sam_align(object):
    def __init__(self):
        self.qname = None
        self.flag = None
        self.rname = None
        self.pos = None
        self.mapq = None
        self.cigar = None
        self.rnext = None
        self.pnext = None
        self.tlen = None
        self.seq = None
        self.qual = None
        self.nm = None

        # secondary alignment position
        self.sa = False
        self.sa_rname = None
        self.sa_pos = None
        self.sa_strand = None
        self.sa_cigar = None
        self.sa_mapq = None
        self.sa_nm = None

    def parse_line(self, text):
        self.text = text.strip()
        sp = self.text.split("\t")

        self.qname = sp[0]
        self.flag = str(bin(int(sp[1])))
        self.rname = sp[2]
        self.pos = int(sp[3])
        self.mapq = int(sp[4])
        self.cigar = sp[5]
        self.rnext = sp[6]
        self.pnext = sp[7]
        self.tlen = int(sp[8])
        self.seq = sp[9]
        self.qual = sp[10]
        #self.nm = int(sp[12].split(':')[2])

        if len(sp) > 16:
            self.sa = True
            sa_sp = sp[-1].split(':')[2].split(',')
            self.sa_rname = sa_sp[0]
            self.sa_pos = int(sa_sp[1])
            self.sa_strand = sa_sp[2]
            self.sa_cigar = sa_sp[3]
            self.sa_mapq = sa_sp[4]
            self.sa_nm = sa_sp[5]


def text_to_sam(text):
    result = sam_align()
    result.parse_line(text)
    return result
