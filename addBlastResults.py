import pandas as pd

xls_blast = pd.ExcelFile('blast_results.xlsx')
xls_intersections = pd.ExcelFile('bedtools_intersections.xlsx')
blast_mirdeep = pd.read_excel(xls_blast, "mirdeep")
intersections_mirdeep = pd.read_excel(xls_intersections, "miRDeep")
"""
example for mirdeep:
for each entry in mirdeep table,
isolate the ID with str.split and str.replace("ID=", ''),
and find the ID in the relevant column in blast results. Use a vectorized method like loc, i did it before.
after locating the row in blast, add to intersect the columns: subject accession, alignment length, e value
save_csv the new intersect table.

after that works, try adding blast to both miRdeep and sRNAbench in every file.
"""

print(intersections_mirdeep.head())
print(blast_mirdeep.head())
for index, row in intersections_mirdeep.iterrows():
    if index == 0:
        continue
    ID = row[8].split(';')[0].replace("ID=", '') + '|'
    print("ID:", ID)
    seed = row[8].split(';')[3]
    print("SEED:", seed)
    mask = blast_mirdeep.loc[blast_mirdeep['query accession'].str.contains(ID)]
    # & blast_mirdeep['query accession'].str.contains(seed)
    print(mask)
    # mask2 = mask.loc[mask['query accession'].str.contains(seed)]
    # print(mask2)


intersections_mirdeep.to_csv('intersections_mirdeep.csv', sep='\t')
