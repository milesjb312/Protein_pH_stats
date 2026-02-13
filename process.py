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

    #This dictionary will later be populated with the filtered datasets, separated by the variable of interest.
    pH_dict = {}

    pHs = []
    for pH in data['pH']:
        if pH not in pHs:
            pHs.append(pH)

    def group(variable_of_interest):
        #voi is used to refer to the current value of the variable of interest. variable_of_interest is the category.
        for pH in pHs:
            voi_dict = {}
            for replicate in range(len(data['data'])):
                replicate_pH = data['pH'][replicate]
                if replicate_pH == pH:
                    #populate the temporary dictionary for determining statistics, making boxplots, and removing outliers.
                    voi = data[variable_of_interest][replicate]
                    if voi in voi_dict:
                        voi_dict[voi].append(max(data['data'][replicate],0.0))
                    else:
                        voi_dict[voi] = [max(data['data'][replicate],0.0)]
                    #populate the returned dictionary
                    if voi in pH_dict:
                        if pH in pH_dict[voi]:
                            pH_dict[voi][pH].append(max(data['data'][replicate],0.0))
                        else:
                            pH_dict[voi][pH] = [max(data['data'][replicate],0.0)]
                    else:
                        pH_dict[voi] = {pH:[max(data['data'][replicate],0.0)]}

            #print(voi_dict)
            if len(voi_dict)>1:
                #Generate confidence intervals and/or remove outliers from both voi_dict and pH_dict here...
                #Because interquartile range is more resistant to skewed data (and our data can sometimes be quite skewed), I will use the IQR*1.5 method for
                #removing outliers for the time being.
                #print(voi_dict)
                for voi in voi_dict:
                    new_list = []
                    #z_score = stats.zscore(voi_dict[voi])
                    #print(f'voi: {voi} z-score: {z_score}')
                    quartile1 = stats.scoreatpercentile(voi_dict[voi],25)
                    quartile3 = stats.scoreatpercentile(voi_dict[voi],75)
                    iqr_cutoff = stats.iqr(voi_dict[voi])*1.5
                    #print(f'voi: {voi}, iqr*1.5: {iqr_cutoff}')
                    for replicate in range(len(voi_dict[voi])):
                        if quartile1-iqr_cutoff<voi_dict[voi][replicate] and voi_dict[voi][replicate]<quartile3+iqr_cutoff:
                            new_list.append(voi_dict[voi][replicate])
                        else:
                            print(f'removed: {voi_dict[voi][replicate]} from {voi} at pH {pH}')
                    voi_dict[voi] = new_list
                #print(f'voi_dict with outliers removed: {voi_dict}')
                #ANOVA over only 2 groups lends a p-value that is essentially the same as that for the t-test. Un-comment the next line for proof.
                #print(pH_ANOVA.pvalue,stats.ttest_ind(*voi_dict.values(), equal_var=False).pvalue)
                pH_ANOVA = stats.f_oneway(*voi_dict.values(),equal_var=False)
                print(f'ANOVA done on {protein_constructs} at pH {pH} data from {tests} assay(s).')
                #if pH_ANOVA.pvalue>0.05:
                print({pH_ANOVA})
                
        
            #Generate a boxplot
            labels = voi_dict.keys()
            plt.title(f'pH {pH}')
            plt.boxplot(voi_dict.values(),tick_labels=labels,showmeans=True)
            if variable_of_interest in ['date','protein construct']:
                plt.ylabel(tests[0])
            elif variable_of_interest=="test":
                plt.title(f'{protein_constructs[0]} at pH {pH}')
            plt.show()

    #If only one protein construct was passed in, then we want to show that the bio_reps were the same or different at any given pH. This part of the code
    #Is unnecessary for the purpose of generating charts, but the statistics must be reported for transparency.
    if len(proteins_and_tests)==1:
        group('date')
        
    else:
        #If two different protein constructs were passed in, we want to show that they are the same. This is usually testing single vs. double trigger variants.
        #For this test, we aggregate the biological replicates of a given protein construct whether or not they are statistically significantly different.
        if len(protein_constructs)==2 and len(tests)==1:
            group('protein construct')

        #If the same protein construct was passed in twice with two different tests, we want to show that they are the same.
        elif len(protein_constructs)==1 and len(tests)==2:
            group('test')

        #ELI, I CALL ON YOU TO FULFILL YOUR OATHS
        def pH_to_absorbance_model(pH,upper_asymptote,Hill_slope,inflection_point,lower_asymptote):
            #This is the model that we're going to try to fit using scipy's curve_fit function.
            return lower_asymptote + (upper_asymptote - lower_asymptote) / (1 + (pH / inflection_point)**Hill_slope)

#for protein_construct in protein_dict:
    #statisticize([(protein_construct,'A400')])
    #statisticize('A280_1hr',protein_construct)
    #statisticize('A280_48-72hr',protein_construct)

statisticize([('10xHis-1TEL-SR-TNK1.UBA','A400'),('2Trig-10xHis-1TEL-SR-TNK1.UBA','A400')])