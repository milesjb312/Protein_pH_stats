import csv
from scipy import stats
from scipy.optimize import curve_fit

#ANOVA must be performed on all groupings of data where there are 3 or more samples, and when it shows that there is a significant difference in means,
#paired T-tests must be performed to see which samples are significantly different.

#All of the data for a given construct at a given pH on the same date must be put into a list.
#Lists must then be compared between proteins of the same construct at the same pH for every different date to see if the date causes a significant change.
#The double and single-trigger variants must also be compared to see if that treatment causes a significant change.

#4-Parameter logistic regression modeling is the closest fit to our data.

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
num_passed = 0

def statisticize(test,protein_construct,protein_construct2=None):
    """
    test is the type of test you performed, referred to by the variable being measured, or the column name in the csv. Ie.: 'A400', 'A280_1hr'
    protein_construct is the dictionary key that refers to the particular protein construct of interest.
    protein_construct2 is an optional second protein construct that you want to compare to the first one.
    So, yes, you do need to run this function once for each type of protein construct and variable combination.
    """
    #Passed data 
    passed_data = {'date':[],'pH':[],test:[]}
    data_and_stat = {'date':[],'pH':[],test:[],'passed':[]}
    if protein_construct2==None:
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
                    if test in ['A280_1hr','A280_48-72hr']:
                        #Divide by the theoretical_max_concentration to normalize between replicates with different stock protein concentrations.
                        #NOTE: this may be inappropriate, but we don't know.
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

                            print(f'ANOVA passed for: {protein_construct} {date} {test} at pH {pH}')
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

        else:
            print(f'There was only one replicate performed for {protein_construct} {test}.')

for protein_construct in protein_dict:
    statisticize('A400',protein_construct)
    #statisticize('A280_1hr',protein_construct)
    #statisticize('A280_48-72hr',protein_construct)

print(num_passed)