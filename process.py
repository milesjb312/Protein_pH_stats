import csv

#ANOVA must be performed on all groupings of data where there are 3 or more samples, and when it shows that there is a significant difference in means,
#paired T-tests must be performed to see which samples are significantly different.

#All of the data for a given construct at a given pH on the same date must be put into a list.
#Lists must then be compared between proteins of the same construct at the same pH for every different date to see if the date causes a significant change.
#The double and single-trigger variants must also be compared to see if that treatment causes a significant change.

protein_dict = {}

with open('data.csv',newline="") as csvfile:
    fieldnames = ['Protein Construct','Date','Stock Concentration mg/mL','pH','A400','A280_1hr','A280_48-72hr']
    reader = csv.DictReader(csvfile,fieldnames=fieldnames)
    next(reader)
    for row in reader:
        if row['Protein Construct'] not in protein_dict:
            protein_dict[row['Protein Construct']] = {}
        if row['Date']+"_"+row['Stock Concentration mg/mL']+'_mg/mL_pH_'+row['pH'] not in protein_dict[row['Protein Construct']]:
            protein_dict[row['Protein Construct']][row['Date']+"_"+row['Stock Concentration mg/mL']+'_mg/mL_pH_'+row['pH']] = {'A400':[],'A280_1hr':[],'A280_48-72hr':[]}
        if row['A400']!="":
            protein_dict[row['Protein Construct']][row['Date']+"_"+row['Stock Concentration mg/mL']+'_mg/mL_pH_'+row['pH']]['A400'].append(float(row['A400']))
        if row['A280_1hr']!="":
            protein_dict[row['Protein Construct']][row['Date']+"_"+row['Stock Concentration mg/mL']+'_mg/mL_pH_'+row['pH']]['A280_1hr'].append(float(row['A280_1hr']))
        if row['A280_48-72hr']!="":
            protein_dict[row['Protein Construct']][row['Date']+"_"+row['Stock Concentration mg/mL']+'_mg/mL_pH_'+row['pH']]['A280_48-72hr'].append(float(row['A280_48-72hr']))

print(protein_dict)