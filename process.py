import csv
from scipy import stats
from scipy.optimize import curve_fit

#First, place all the data within a dictionary, separating by protein construct, then by date/concentration/pH, then by test-type.
protein_dict = {}

with open('data.csv',newline="") as csvfile:
    fieldnames = ['Protein Construct','Date','Stock Concentration mg/mL','pH','A400','A280_1hr','A280_48-72hr']
    reader = csv.DictReader(csvfile,fieldnames=fieldnames)
    next(reader)
    for row in reader:
        protein_construct = row['Protein Construct'].strip()
        if protein_construct not in protein_dict:
            print(f'{protein_construct} placed in dictionary.')
            protein_dict[protein_construct] = {}
        if row['Date']+"_"+row['Stock Concentration mg/mL']+'_mg/mL_pH_'+row['pH'] not in protein_dict[protein_construct]:
            protein_dict[protein_construct][row['Date']+"_"+row['Stock Concentration mg/mL']+'_mg/mL_pH_'+row['pH']] = {'A400':[],'A280_1hr':[],'A280_48-72hr':[]}
        if row['A400']!="":
            protein_dict[protein_construct][row['Date']+"_"+row['Stock Concentration mg/mL']+'_mg/mL_pH_'+row['pH']]['A400'].append(float(row['A400']))
        if row['A280_1hr']!="":
            protein_dict[protein_construct][row['Date']+"_"+row['Stock Concentration mg/mL']+'_mg/mL_pH_'+row['pH']]['A280_1hr'].append(float(row['A280_1hr']))
        if row['A280_48-72hr']!="":
            protein_dict[protein_construct][row['Date']+"_"+row['Stock Concentration mg/mL']+'_mg/mL_pH_'+row['pH']]['A280_48-72hr'].append(float(row['A280_48-72hr']))

#print(protein_dict)
num_passed = 0

def statisticize(tests:list,protein_constructs:list):
    """
    Normally, you should pass only one test and one protein construct into tests and protein constructs, in the form of a list.
    This shows which sets of data from different dates but the same pH and test are significantly different.
    'tests' are the types of tests you performed, referred to by the variable being measured, or the column name in the csv. Ie.: 'A400', 'A280_1hr'
    'protein_constructs' are the protein_dict dictionary keys that refer to the particular protein constructs of interest.

    So, yes, you do need to run this function once for each type of protein construct and variable combination.
    """
    for test in tests:
        for protein_construct in protein_constructs:
            #Passed data 
            passed_data = {'date':[],'pH':[],test:[]}
            data_and_stat = {'date':[],'pH':[],test:[],'passed':[]}
            
            data = {'date':[],'pH':[],test:[]}
            dates = []
            pHs = []
            for bio_rep in protein_dict[protein_construct]:
                if len(protein_dict[protein_construct][bio_rep][test])!=0:
                    #if the biological replicate actually had this test done, aggregate its data into the temporary dataframe replicate by replicate, conserving their order.
                    for replicate in protein_dict[protein_construct][bio_rep][test]:

                        date = bio_rep.split("_")[0]
                        if date not in data['date']:
                            dates.append(date)
                        data['date'].append(bio_rep.split("_")[0])

                        pH = bio_rep.split("_")[-1]
                        if pH not in data['pH']:
                            pHs.append(pH)
                        data['pH'].append(pH)

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
                            data[test].append((theoretical_max_concentration-replicate)/theoretical_max_concentration)
                        else:
                            data[test].append(replicate/theoretical_max_concentration)
            #do ANOVA on each set of pH values.
            if len(dates)>1:
                for pH in pHs:
                    #For each pH, make a temporary dictionary that splits the replicates by biological replicate.
                    pH_date_dict = {}
                    for replicate in range(len(data['pH'])):
                        if data['pH'][replicate] == pH:
                            date = data['date'][replicate]
                            if date in pH_date_dict:
                                pH_date_dict[date].append(data[test][replicate])
                            else:
                                pH_date_dict[date] = [data[test][replicate]]
                    #print(pH_date_dict)
                    pH_ANOVA = stats.f_oneway(*pH_date_dict.values(),equal_var=False)
                    print(f'{protein_construct} {pH} {test}')
                    #Generate standard error here...
                    if len(dates)>2:
                        #ANOVA over only 2 groups lends a p-value that is essentially the same as that for the t-test. Un-comment the next line for proof.
                        #print(pH_ANOVA.pvalue,stats.ttest_ind(*pH_date_dict.values(), equal_var=False).pvalue)
                        if pH_ANOVA.pvalue>=0.05:
                            for date in pH_date_dict:
                                for replicate in pH_date_dict[date]:
                                    passed_data['pH'].append(pH)
                                    passed_data[test].append(replicate)

                                    data_and_stat['date'].append(date)
                                    data_and_stat['pH'].append(pH)
                                    data_and_stat[test].append(replicate)
                                    data_and_stat['passed'].append(True)

                            print(f'ANOVA passed for: {protein_construct} {test} at pH {pH}')
                            global num_passed
                            num_passed+=1

                        else:
                            data_and_stat['date'].append(date)
                            data_and_stat['pH'].append(pH)
                            data_and_stat[test].append(replicate)
                            data_and_stat['passed'].append(False)

                #return passed_data, data_and_stat OR create the charts needed within this function itself (which is what I would suggest)
                #ELI, I CALL ON YOU TO FULFILL YOUR OATHS
                #Do we want to stick to a logistic model?
                #A logistic model claims that the polymer does not have variable stable lengths at different pHs, but rather it forms
                #under all conditions below a certain pH and does not form above that pH.
                #Alternatively, we can do a polynomial fit, which can produce several inflection points.
                #For now, I will stick to a 4-parameter logistic regression.

                def pH_to_absorbance_model(pH,upper_asymptote,Hill_slope,inflection_point,lower_asymptote):
                    #This is the model that we're going to try to fit using scipy's curve_fit function.
                    return lower_asymptote + (upper_asymptote - lower_asymptote) / (1 + (pH / inflection_point)**Hill_slope)

                #Here we would pass the pH_to_absorbance_model along with our passed_data into a curve_fit. This is problematic, since we only have 16
                #sets of data that pass the ANOVA (meaning all biological replicates at a particular pH point had the same mean) and they come from only
                #5 of the constructs. They also don't surround the theoretical inflection point in any but a single construct. So I won't be doing this.

                # We have two other options. One: pass all the data.
                #We generate curves for each biological replicate. This would probably better show the differences
                #between them. We could even then show the average curve as the mean of all the curves for a particular construct. This is the most elegant
                #and transparent approach, in my opinion. The mean curve would just take the 4-parameter logistic variables and do the mean of each of them.

                #Alternatively, we could try something a little harder.
                #We group the pH points that are below and above a theoretical inflection point, testing all possible points for each construct.
                #Then, we do ANOVA on the combined groups. This could let us, for instance, combine all the data from pH points 4, 4.5, and 5 and assign them
                # to the average pH point of the group. This generates a wide spread of readings before the inflection point that are less likely to be 
                #statistically different from the readings from other dates by sacrificing precision (the model doesn't know the high and low pH points precisely).
                #Any sets of two groups that pass the ANOVA could theoretically produce the closest curve. All sets that pass would be used to generate a curve.
                #The code will be a lot harder, and I'm not convinced it would perform better than the previous option.

            else:
                print(f'There was only one replicate performed for {protein_construct} {test}.')

for protein_construct in protein_dict:
    statisticize('A400',protein_construct)
    #statisticize('A280_1hr',protein_construct)
    #statisticize('A280_48-72hr',protein_construct)

print(num_passed)