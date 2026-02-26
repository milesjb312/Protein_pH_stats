import csv
from scipy import stats
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from matplotlib.legend_handler import HandlerTuple
import numpy as np

#1)Pull specified constructs and their corresponding specified data from csv and put into protein_dict by calling the master function 'statisticize'. Ex.:
    #statisticize([('Construct1','A400'),('Construct2','A400')])
#2)Repeat steps 3-4 for each construct asked for.
#3)Either combine or don't all of the biological replicates (combine bool==True or False). Populate the grouped_dict with this.
#4)Determine medians, try to fit a curve, and plot.

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
            if '*' in row['A400']:
                protein_dict[protein_construct][row['Date']+"_"+row['Stock Concentration mg/mL']+'_mg/mL_pH_'+row['pH']]['A400'].append((row['A400']))
            else:
                protein_dict[protein_construct][row['Date']+"_"+row['Stock Concentration mg/mL']+'_mg/mL_pH_'+row['pH']]['A400'].append(float(row['A400']))
        if row['A280_1hr']!="":
            if '*' in row['A280_1hr']:
                protein_dict[protein_construct][row['Date']+"_"+row['Stock Concentration mg/mL']+'_mg/mL_pH_'+row['pH']]['A280_1hr'].append((row['A280_1hr']))
            else:
                protein_dict[protein_construct][row['Date']+"_"+row['Stock Concentration mg/mL']+'_mg/mL_pH_'+row['pH']]['A280_1hr'].append(float(row['A280_1hr']))
        if row['A280_48-72hr']!="":
            if '*' in row['A280_48-72hr']:
                protein_dict[protein_construct][row['Date']+"_"+row['Stock Concentration mg/mL']+'_mg/mL_pH_'+row['pH']]['A280_48-72hr'].append((row['A280_48-72hr']))
            else:
                protein_dict[protein_construct][row['Date']+"_"+row['Stock Concentration mg/mL']+'_mg/mL_pH_'+row['pH']]['A280_48-72hr'].append(float(row['A280_48-72hr']))

def statisticize(proteins_and_tests:list,boxplot=False,plot=False,combine=False,require_matched_groups_var=False):
    """
    Normally, you should pass only one test and one protein construct into 'proteins_and_tests', in the form of a list of tuples, where the
    tuple looks like: (protein_construct,test)
    'protein_constructs' are the protein_dict dictionary keys that refer to the particular protein constructs of interest.
    'tests' are the types of tests you performed, referred to by the variable being measured, or the column name in the csv. Ie.: 'A400', 'A280_1hr'    
    """
    #The first half of this function is just a filter pulling applicable data from the protein_dict dictionary.
    data = {'protein construct':[],'date':[],'pH':[],'test':[],'data':[],'tailored_data':[]}
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

                    if type(replicate)==float:
                        if test in ['A280_1hr','A280_48-72hr']:
                            #For precipitation assays, find the amount of protein left as a percentage.
                            data['tailored_data'].append((theoretical_max_concentration-replicate)/theoretical_max_concentration)
                            data['data'].append((theoretical_max_concentration-replicate)/theoretical_max_concentration)
                        else:
                            data['tailored_data'].append(replicate/theoretical_max_concentration)
                            data['data'].append(replicate/theoretical_max_concentration)
                    else:
                        if test in ['A280_1hr','A280_48-72hr']:
                            #For precipitation assays, find the amount of protein left as a percentage.
                            data['data'].append((theoretical_max_concentration-float(replicate.strip("*")))/theoretical_max_concentration)
                            data['tailored_data'].append('')
                        else:
                            data['data'].append(float(replicate.strip("*"))/theoretical_max_concentration)
                            data['tailored_data'].append('')

    #--------------------Now that the main filtering is done, it's time to do ANOVA or Kruskal-Wallis on each set of pH values--------------------------------#
    #The first step is to further filter the data, moving it into dictionaries that denote the groups of interest. This is dependent on the number of
    #protein constructs and tests passed into this function.

    #This dictionary will later be populated with the filtered datasets, separated by date or by a combined grouping of some variable of interest.
    grouped_dict = {}

    pHs = []
    for pH in data['pH']:
        if pH not in pHs:
            pHs.append(pH)

    def group(combine):
        #If you want to combine the replicates:
        if combine:
            #If only one protein construct was passed in, then we want to show that the bio_reps were the same or different at any given pH. This part of the code
            #Is unnecessary for the purpose of generating charts, but the statistics must be reported for transparency.
            if len(proteins_and_tests)==1:
                variable_of_interest = 'date'    
            else:
                #If two different protein constructs were passed in, we want to show that they are the same. This is usually testing single vs. double trigger variants.
                #For this test, we aggregate the biological replicates of a given protein construct whether or not they are statistically significantly different.
                if len(protein_constructs)==2 and len(tests)==1:
                    variable_of_interest = 'protein construct'

                #If the same protein construct was passed in twice with two different tests, we want to show that they are the same.
                elif len(protein_constructs)==1 and len(tests)==2:
                    variable_of_interest = 'test'

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
                        if voi in grouped_dict:
                            if pH in grouped_dict[voi]:
                                grouped_dict[voi][pH].append(max(data['data'][replicate],0.0))
                            else:
                                grouped_dict[voi][pH] = [max(data['data'][replicate],0.0)]
                        else:
                            grouped_dict[voi] = {pH:[max(data['data'][replicate],0.0)]}

                #Generate confidence intervals here...
                pH_date_kruskal = stats.kruskal(*voi_dict.values())
                #print(f'{pH}: {pH_date_kruskal}')
                #pH_date_ANOVA = stats.f_oneway(*voi_dict.values(),equal_var=False)
                #print(pH_date_ANOVA)
                #ANOVA over only 2 groups lends a p-value that is essentially the same as that for the t-test. Un-comment the next line for proof.
                #print(pH_ANOVA.pvalue,stats.ttest_ind(*voi_dict.values(), equal_var=False).pvalue)
                #print(f'ANOVA done on {protein_constructs} at pH {pH} data from {tests} assay(s).')
                #if pH_ANOVA.pvalue>0.05:
                #print({pH_ANOVA})

                #Generate a boxplot
                if boxplot==True:
                    labels = voi_dict.keys()
                    plot = plt.boxplot(voi_dict.values(),tick_labels=labels,showmeans=True)
                    if len(proteins_and_tests)==1:
                        plt.title(f'{protein_constructs[0]} {tests[0]} at pH {pH}')
                        plt.ylabel(tests[0])
                    else:
                        #If two different protein constructs were passed in, we want to show that they are the same. This is usually testing single vs. double trigger variants.
                        #For this test, we aggregate the biological replicates of a given protein construct whether or not they are statistically significantly different.
                        if len(protein_constructs)==2 and len(tests)==1:
                            plt.title(f'{protein_constructs[0]} vs. {protein_constructs[1]} {tests[0]} at pH {pH}')
                            plt.ylabel(tests[0])
                        #If the same protein construct was passed in twice with two different tests, we want to show that they are the same.
                        elif len(protein_constructs)==1 and len(tests)==2:
                            plt.title(f'{protein_constructs[0]} {tests[0]} vs. {tests[1]} at pH {pH}')
                            plt.ylabel(tests[0])
                    plt.show()

        #If you don't want to combine the replicates:
        else:
            for pH in pHs:
                for replicate in range(len(data['tailored_data'])):
                    if type(data['tailored_data'][replicate])==float:
                        if pH == data['pH'][replicate]:
                            if len(tests)==2:
                                construct_and_date = data['protein construct'][replicate]+"_"+data['date'][replicate]+"_"+data['test'][replicate]
                            else:
                                construct_and_date = data['protein construct'][replicate]+"_"+data['date'][replicate]
                            if construct_and_date in grouped_dict:
                                if pH in grouped_dict[construct_and_date]:
                                    grouped_dict[construct_and_date][pH].append(max(data['tailored_data'][replicate],0.0))
                                else:
                                    grouped_dict[construct_and_date][pH] = [max(data['tailored_data'][replicate],0.0)]
                            else:
                                grouped_dict[construct_and_date] = {pH:[max(data['tailored_data'][replicate],0.0)]}

    def require_matched_groups():
        """This function is meant to be used when comparing two tests of the same biological replicate. It gets rid of biological replicates with only one of
        the tests."""
        new_grouped_dict = {}
        date_counts = {}
        for construct_and_date in grouped_dict:
            date=construct_and_date.split("_")[1]
            if date in date_counts:
                date_counts[date]+=1
            else:
                date_counts[date]=1
        for construct_and_date in grouped_dict:
            if date_counts[construct_and_date.split("_")[1]]>=2:
                new_grouped_dict[construct_and_date] = grouped_dict[construct_and_date]
        return new_grouped_dict

    group(combine)

    if require_matched_groups_var:
        grouped_dict = require_matched_groups()

    def pH_to_absorbance_model_4pl(pH,upper_asymptote,Hill_slope,inflection_point,lower_asymptote):
        #This is the model that we're going to try to fit using scipy's curve_fit function.
        return lower_asymptote + (upper_asymptote - lower_asymptote) / (1 + (pH / inflection_point)**Hill_slope)

    if plot==True:
        handles = []
        labels = []
        inflection_handle = None
        inflection_label = "Inflection Point"
        keys_1trig = []
        keys_2trig = []
        for protein in grouped_dict:
            if "2Trig" in (protein.split("_")[0]):
                    keys_2trig.append(protein)
            else:
                keys_1trig.append(protein)
        def shade_list(cmap, n, lo=0.25, hi=0.9):
            if n<= 1:
                return [cmap((lo + hi) / 2)]
            return [cmap(x) for x in np.linspace(lo, hi, n)]
        blue_shades = shade_list(plt.cm.Blues, len(keys_1trig))
        orange_shades = shade_list(plt.cm.Oranges, len(keys_2trig))
        color_map = {}
        for protein, color in zip(sorted(keys_1trig), blue_shades):
            color_map[protein] = color
        for protein, color in zip(sorted(keys_2trig), orange_shades):
            color_map[protein] = color
        inflection_points_1trig = []
        inflection_points_2trig = []
        for construct_and_date in grouped_dict:
            #Somewhere below this line, we would incorporate the Confidence Interval Band.
            #What would realistically be a 95% confidence interval here? We don't know the shape of the population, so we could assume a normal distribution,
            #but it would be somewhat meaningless. If we did so, the confidence interval would be constructed using the standard error of either each
            #biological replicate at a single pH or the combined biological replicates for a single construct at a single pH. Because the charts would be too
            #messy with a CI band and 6 or more curves, I propose that we should only include them in the combined charts, and only if we feel it is entirely
            #necessary.
            """
            plt.fill_between(
                fitted_curve,
                fitted_curve - lower_confidence_limit,
                fitted_curve + upper_confidence_limit,
                color=color_band,
                alpha=0.45,
                label="95% CI band",
            )"""
            #HERE IS WHERE TO WORK NEXT! All of the plots comparing single to double trigger show both the curves and scatter dots. We need to now make it so that anything
            #combined doesn't include the scatter plot dots.
            try:
                pH_values = []
                medians = []
                for pH in grouped_dict[construct_and_date]:
                    if len(grouped_dict[construct_and_date][pH])!=0:
                        pH_values.append(float(pH))
                    #mean = sum(grouped_dict[construct_and_date][pH])/len(grouped_dict[construct_and_date][pH])
                    #means.append(mean)
                    median = np.median(grouped_dict[construct_and_date][pH])
                    medians.append(median)
                pH_linspace = np.linspace(min(pH_values),max(pH_values),400)
                #print(f'medians: {medians}')
                #stats.mannwhitneyu()
                initial_guesses = [max(medians),-3,5.5,min(medians)]
                best_fit_parameters,covariance_matrix = curve_fit(pH_to_absorbance_model_4pl,pH_values,medians,p0=initial_guesses)
                #print(f'best fit parameters for {construct_and_date}: {best_fit_parameters}')
                if "/" in construct_and_date:
                    if "2Trig" in construct_and_date and best_fit_parameters[2] > 4.0:
                        inflection_points_2trig.append(best_fit_parameters[2])
                    elif "2Trig" in construct_and_date and best_fit_parameters[2] < 4.0:
                        inflection_points_2trig.append(float(4.0))
                    elif not "2Trig" in construct_and_date and best_fit_parameters[2] < 4.0:
                        inflection_points_1trig.append(float(4.0))
                    else:
                        inflection_points_1trig.append(best_fit_parameters[2])
                #if "/" not in construct_and_date:

                #don't plot the curve_fit if the inflection point is below 4.0
                if best_fit_parameters[2] < 4.0:
                    raise ValueError(f"Calculated inflection point ({round(best_fit_parameters[2],2)}) for {construct_and_date} below 4.0")
                #print(best_fit_parameters)
                if "2Trig" in construct_and_date:
                    plot = plt.plot(pH_linspace,pH_to_absorbance_model_4pl(pH_linspace,*best_fit_parameters),linestyle="--", color=color_map[construct_and_date])
                else:
                    plot = plt.plot(pH_linspace,pH_to_absorbance_model_4pl(pH_linspace,*best_fit_parameters), color=color_map[construct_and_date])
                line_handle = plot[0]
                construct_pHs = []
                construct_data = []
                for pH in grouped_dict[construct_and_date]:
                    for replicate in grouped_dict[construct_and_date][pH]:
                        construct_pHs.append(float(pH))
                        construct_data.append(replicate)
                scatter_handle = None
                if "/" in construct_and_date:
                    if "2Trig" in construct_and_date:
                        scatter_handle = plt.scatter(construct_pHs,construct_data,marker="^", color=color_map[construct_and_date])
                    else:
                        scatter_handle = plt.scatter(construct_pHs,construct_data, color=color_map[construct_and_date])
                else:
                    #The asterisk unpacks all of the variables in bets_fit_parameters so that the model has all the constants it needs to calculate the y value 
                    #at the inflection point.
                    inflection_point_marker = plt.scatter(best_fit_parameters[2], pH_to_absorbance_model_4pl(best_fit_parameters[2],*best_fit_parameters), color='Black', marker='x', s=100, 
                                linewidths=2.5, label="Inflection Point", zorder=10)
                    if inflection_handle is None:
                        inflection_handle = inflection_point_marker
                #Legend Entry - if the scatter_handle is something then we use two symbols in the legend but if it is None then we only use one symbol.
                if scatter_handle is not None:
                    handles.append((line_handle, scatter_handle))
                else:
                    handles.append(line_handle)
                labels.append(construct_and_date)
            except Exception as e:
                print(f'The curve fit for {construct_and_date} failed due to: {e}. Plotting points instead of fitted curve...')
                construct_pHs = []
                construct_data = []
                for pH in grouped_dict[construct_and_date]:
                    for replicate in grouped_dict[construct_and_date][pH]:
                        construct_pHs.append(float(pH))
                        construct_data.append(replicate)
                if "2Trig" in construct_and_date:
                    plot = plt.scatter(construct_pHs,construct_data,marker="^", color=color_map[construct_and_date])
                else:
                    plot = plt.scatter(construct_pHs,construct_data, color=color_map[construct_and_date])
                handles.append(plot)
                labels.append(construct_and_date)

        if len(proteins_and_tests)==1:
            plt.title(f'{protein_constructs[0]} {tests[0]}')
            plt.ylabel(tests[0])
            plt.xlabel('pH')
        else:
            #If two different protein constructs were passed in, we want to show that they are the same. This is usually testing single vs. double trigger variants.
            #For this test, we aggregate the biological replicates of a given protein construct whether or not they are statistically significantly different.
            if len(protein_constructs)==2 and len(tests)==1:
                plt.title(f'{protein_constructs[0]} vs. {protein_constructs[1]} {tests[0]}')
                plt.ylabel(tests[0])
                plt.xlabel('pH')
            #If the same protein construct was passed in twice with two different tests, we want to show that they are the same.
            elif len(protein_constructs)==1 and len(tests)==2:
                plt.title(f'{protein_constructs[0]} {tests[0]} vs. {tests[1]}')
                plt.ylabel(tests[0])
                plt.xlabel('pH')
        legend_handles = handles.copy()
        legend_labels = labels.copy()
        if inflection_handle is not None:
            legend_handles.append(inflection_handle)
            legend_labels.append(inflection_label)
        plt.legend(legend_handles,legend_labels,handler_map={tuple: HandlerTuple(ndivide=2, pad=1)})
        print(f'Single Trigger Inflection Points: {inflection_points_1trig}')
        print(f'Double Trigger Inflection Points: {inflection_points_2trig}')
        plt.show()

for protein_construct in protein_dict:
    if '2Trig' not in protein_construct and 'Gravity' not in protein_construct:
        #Single vs. Double-trigger A400
        statisticize([(f'{protein_construct}','A400'),(f'2Trig-{protein_construct}','A400')],boxplot=False,plot=True,combine=False)
        statisticize([(f'{protein_construct}','A400'),(f'2Trig-{protein_construct}','A400')],boxplot=False,plot=True,combine=True)

#for protein_construct in protein_dict:
 #   if '2Trig' not in protein_construct and 'Gravity' not in protein_construct:
        #Single vs. Double-trigger A280-1hr
  #      statisticize([(f'{protein_construct}','A280_1hr'),(f'2Trig-{protein_construct}','A280_1hr')],boxplot=False,plot=True,combine=False)
    #    statisticize([(f'{protein_construct}','A280_1hr'),(f'2Trig-{protein_construct}','A280_1hr')],boxplot=False,plot=True,combine=True)

#for protein_construct in protein_dict:
    #if '2Trig' not in protein_construct and 'Gravity' not in protein_construct:
        #Single vs. Double-trigger A280_48-72hr
        #statisticize([(f'{protein_construct}','A280_48-72hr'),(f'2Trig-{protein_construct}','A280_48-72hr')],boxplot=False,plot=True,combine=False)
        #statisticize([(f'{protein_construct}','A280_48-72hr'),(f'2Trig-{protein_construct}','A280_48-72hr')],boxplot=False,plot=True,combine=True)

#for protein_construct in protein_dict:
 #   if '2Trig' not in protein_construct and 'Gravity' not in protein_construct:
  #      #Single trigger A280_1hr vs A280_48-72hr
   #     statisticize([(f'{protein_construct}','A280_1hr'),(f'{protein_construct}','A280_48-72hr')],boxplot=False,plot=True,combine=False,require_matched_groups_var=True)

#Single trigger TV-vWA A400
#statisticize([('10xHis-1TEL-TV-vWA','A400'),('10xHis-1TEL-TV-vWA (Gravity)','A400')],plot=True,combine=False)
#statisticize([('10xHis-1TEL-TV-vWA','A400'),('10xHis-1TEL-TV-vWA (Gravity)','A400')],plot=True,combine=True)
