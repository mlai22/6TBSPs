import swalign
# choose your own values here… 2 and -1 are common.
match = 4
mismatch = -1
scoring = swalign.NucleotideScoringMatrix(match, mismatch)

sw = swalign.LocalAlignment(scoring, gap_penalty=-3, gap_extension_penalty=-1)  # you can also choose gap penalties, etc...
# alignment = sw.align('ACACACTA','AGCACACA')
# alignment = sw.align('MISLIAALAVDRVIGMENAMPFNLPADLAWFKRNTLDKPVIMGRHTWESIG','SLNCIVAVSQNMGIGKNGDLPWPPLRNEFRYFQRMTTTSSVEGKQNLVIMGKKTWFSIPE')
alignment = sw.align('ATTGAGC', 'ATGC')
alignment.dump()


if __name__ == "__main__":
    pass