import csv
from scipy import stats
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

#First, place all the data within a dictionary, separating by protein construct, then by date/concentration/pH, then by test-type.
protein_dict = {}

with open('data.csv',newline="") as csvfile:
    fieldnames = ['Protein Construct','Date','Stock Concentration mg/mL','pH','A400','A280_1hr','A280_48-72hr']
    reader = csv.DictReader(csvfile,fieldnames=fieldnames)
    next(reader)
    for row in reader:
        protein_construct = row['Protein Construct'].strip()
        if protein_construct not in protein_dict:
            #print(f'{protein_construct} placed in dictionary.')
            protein_dict[protein_construct] = {}
        if row['Date']+"_"+row['Stock Concentration mg/mL']+'_mg/mL_pH_'+row['pH'] not in protein_dict[protein_construct]:
            protein_dict[protein_construct][row['Date']+"_"+row['Stock Concentration mg/mL']+'_mg/mL_pH_'+row['pH']] = {'A400':[],'A280_1hr':[],'A280_48-72hr':[]}
        if row['A400']!="":
            protein_dict[protein_construct][row['Date']+"_"+row['Stock Concentration mg/mL']+'_mg/mL_pH_'+row['pH']]['A400'].append(float(row['A400']))
        if row['A280_1hr']!="":
            protein_dict[protein_construct][row['Date']+"_"+row['Stock Concentration mg/mL']+'_mg/mL_pH_'+row['pH']]['A280_1hr'].append(float(row['A280_1hr']))
        if row['A280_48-72hr']!="":
            protein_dict[protein_construct][row['Date']+"_"+row['Stock Concentration mg/mL']+'_mg/mL_pH_'+row['pH']]['A280_48-72hr'].append(float(row['A280_48-72hr']))

def statisticize(proteins_and_tests:list):
    """
    Normally, you should pass only one test and one protein construct into 'proteins_and_tests', in the form of a list of tuples, where the
    tuple looks like: (protein_construct,test)
    'protein_constructs' are the protein_dict dictionary keys that refer to the particular protein constructs of interest.
    'tests' are the types of tests you performed, referred to by the variable being measured, or the column name in the csv. Ie.: 'A400', 'A280_1hr'    
    """
    #The first half of this function is just a filter pulling applicable data from the protein_dict dictionary.
    data = {'protein construct':[],'date':[],'pH':[],'test':[],'data':[]}
    protein_constructs = []
    tests = []
    for protein_construct,test in proteins_and_tests:
        if protein_construct not in protein_constructs:
            protein_constructs.append(protein_construct)
        if test not in tests:
            tests.append(test)
        for bio_rep in protein_dict[protein_construct]:
            if len(protein_dict[protein_construct][bio_rep][test])!=0:
                #if the biological replicate actually had this test done, aggregate its data into the dataframe replicate by replicate, conserving their order.
                for replicate in protein_dict[protein_construct][bio_rep][test]:
                    data['protein construct'].append(protein_construct)
                    date = bio_rep.split("_")[0]
                    data['date'].append(date)
                    pH = bio_rep.split("_")[-1]
                    data['pH'].append(pH)
                    data['test'].append(test)

                    concentration = float(bio_rep.split("_")[1])
                    #We dilute our protein by putting 200 uL into 800 uL.
                    theoretical_max_concentration = concentration/5
                    #Divide by the theoretical_max_concentration to normalize between replicates with different stock protein concentrations.
                    #This means that for precipitation assays, we are reporting the percent of remaining soluble protein,
                    #and for turbidity assays, we are reporting A400 absorbance per mg/mL of protein.
                    #Beer's law says that:
                    #Absorbance=extinction_coefficient*column_length*concentration.
                    #This means it should be appropriate for us to divide absorbance by concentration. If we do so, we are effectively reporting on the
                    #extinction coefficient of the precipitated protein. Another word for this is the molar absorptivity. We propose that the higher a
                    #molar absorptivity is, the denser the precipitate is. Essentially, the solubility and turbidity assays are testing entirely
                    #different things.

                    if test in ['A280_1hr','A280_48-72hr']:
                        #For precipitation assays, find the amount of protein left as a percentage.
                        data['data'].append((theoretical_max_concentration-replicate)/theoretical_max_concentration)
                    else:
                        data['data'].append(replicate/theoretical_max_concentration)
            else:
                print(f'{protein_construct} on {bio_rep} has no {test} test data.')
    #--------------------Now that the main filtering is done, it's time to do ANOVA on each set of pH values-------------------------------------------#
    #The first step is to further filter the data, moving it into dictionaries that denote the groups of interest. This is dependent on the number of
    #protein constructs and tests passed into this function.
    #This could probably be done initially, but I can't figure it out.
    pHs = []
    for pH in data['pH']:
        if pH not in pHs:
            pHs.append(pH)

    #If only one protein construct was passed in, then we want to show that the bio_reps were the same or different at any given pH. This part of the code
    #Is unnecessary for the purpose of generating charts, but the statistics must be reported for transparency.
    if len(proteins_and_tests)==1:
        for pH in pHs:
            date_dict = {'date':[]}
            for replicate in range(len(data['data'])):
                replicate_pH = data['pH'][replicate]
                if replicate_pH == pH:
                    date = data['date'][replicate]
                    if date in date_dict:
                        date_dict[date].append(max(data['data'][replicate],0.0))
                    else:
                        date_dict[date] = [max(data['data'][replicate],0.0)]
            #print(date_dict)
            if len(date_dict)>1:
                #ANOVA over only 2 groups lends a p-value that is essentially the same as that for the t-test. Un-comment the next line for proof.
                #print(pH_ANOVA.pvalue,stats.ttest_ind(*pH_date_dict.values(), equal_var=False).pvalue)
                pH_ANOVA = stats.f_oneway(*date_dict.values(),equal_var=False)
                print(f'{protein_constructs[0]} {pH} {test}')
                #if pH_ANOVA.pvalue>0.05:
                print({pH_ANOVA})
                #Generate standard error here...
    
    else:
        #If two different protein constructs were passed in, we want to show that they are the same. This is usually testing single vs. double trigger variants.
        #For this test, we aggregate the biological replicates of a given protein construct whether or not they are statistically significantly different.
        if len(protein_constructs)==2 and len(tests)==1:
            for pH in pHs:
                construct_dict = {}
                for replicate in range(len(data['data'])):
                    replicate_pH = data['pH'][replicate]
                    if replicate_pH == pH:
                        protein_construct = data['protein construct'][replicate]
                        if protein_construct in construct_dict:
                            construct_dict[protein_construct].append(max(data['data'][replicate],0.0))
                        else:
                            construct_dict[protein_construct] = [max(data['data'][replicate],0.0)]

                if len(construct_dict)>1:
                    #ANOVA over only 2 groups lends a p-value that is essentially the same as that for the t-test. Un-comment the next line for proof.
                    #print(pH_ANOVA.pvalue,stats.ttest_ind(*pH_date_dict.values(), equal_var=False).pvalue)
                    pH_ANOVA = stats.f_oneway(*construct_dict.values(),equal_var=False)
                    labels = construct_dict.keys()
                    plt.boxplot(construct_dict.values(),tick_labels=labels)
                    plt.show()
                    #print(f'{protein_constructs} {pH} {tests[0]}')
                    #if pH_ANOVA.pvalue>0.05:
                    print({pH_ANOVA})
                    #Generate standard error here...

        #If the same protein construct was passed in twice with two different tests, we want to show that they are the same.
        elif len(protein_constructs)==1 and len(tests)==2:
            pass

        #ELI, I CALL ON YOU TO FULFILL YOUR OATHS
        def pH_to_absorbance_model(pH,upper_asymptote,Hill_slope,inflection_point,lower_asymptote):
            #This is the model that we're going to try to fit using scipy's curve_fit function.
            return lower_asymptote + (upper_asymptote - lower_asymptote) / (1 + (pH / inflection_point)**Hill_slope)

#for protein_construct in protein_dict:
    #statisticize([(protein_construct,'A400')])
    #statisticize('A280_1hr',protein_construct)
    #statisticize('A280_48-72hr',protein_construct)

statisticize([('1TEL-GG-TNK1.UBA','A400'),('2Trig-1TEL-GG-TNK1.UBA','A400')])