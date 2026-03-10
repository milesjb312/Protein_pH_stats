import csv
from scipy import stats
from scipy.optimize import curve_fit
from scipy.stats import t
import matplotlib.pyplot as plt
from matplotlib.legend_handler import HandlerTuple
import numpy as np

#1)Pull specified constructs and their corresponding specified data from csv and put into protein_dict by calling the master function 'statisticize'. Ex.:
    #statisticize([('Construct1','A400'),('Construct2','A400')])
#2)Repeat steps 3-4 for each construct asked for.
#3)For the purposes of data analysis, only combine biological replicates right before making plots.
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

def statisticize(proteins_and_tests:list,boxplot=False,plot_singly=False,plot=False,require_matched_groups_var=False,tailored=False):
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

    pHs = []
    for pH in data['pH']:
        if pH not in pHs:
            pHs.append(pH)

    if tailored == True:
        data_quality = 'tailored_data'
    else:
        data_quality = 'data'

    # Eli's Notes: 
     # data[data_quality] is the dictionary that contains all of our data
     # len(data[data_quality]) returns the numbers of rows in the dataset
     # range(len(data[data_quality])) creates a sequence of integers starting from 0 and going for as long as len(data[data_quality])
     # i corresponds to one row in the dataset; each row was determined with the range(...) function; eg i will go from 0 to len(data[data_quality])
        # In summary, this for loop iterates through each data point measured on the nanodrop (all of the data points are in the data dictionary) and adds each data point
        # to the new constructs_by_pH_by_date dictionary attached to three keys (construct, pH, date).
    constructs_by_pH_by_date = {}
    for i in range(len(data[data_quality])):
        construct = data['protein construct'][i]
        pH = data['pH'][i]
        date = data['date'][i]
        measurement = data[data_quality][i]
        if construct not in constructs_by_pH_by_date:
            constructs_by_pH_by_date[construct] = {}
        if pH not in constructs_by_pH_by_date[construct]:
            constructs_by_pH_by_date[construct][pH] = {}
        if date not in constructs_by_pH_by_date[construct][pH]:
            constructs_by_pH_by_date[construct][pH][date] = []
        constructs_by_pH_by_date[construct][pH][date].append(measurement)

    """
    for pH in pHs:
        for replicate in range(len(data[data_quality])):
            if type(data[data_quality][replicate])==float:
                if pH == data['pH'][replicate]:
                    if len(tests)==2:
                        construct_and_date = data['protein construct'][replicate]+"_"+data['date'][replicate]+"_"+data['test'][replicate]
                    else:
                        construct_and_date = data['protein construct'][replicate]+"_"+data['date'][replicate]
                    if construct_and_date in grouped_dict:
                        if pH in grouped_dict[construct_and_date]:
                            grouped_dict[construct_and_date][pH].append(max(data[data_quality][replicate],0.0))
                        else:
                            grouped_dict[construct_and_date][pH] = [max(data[data_quality][replicate],0.0)]
                    else:
                        grouped_dict[construct_and_date] = {pH:[max(data[data_quality][replicate],0.0)]}
    """

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

    if require_matched_groups_var:
        grouped_dict = require_matched_groups()

    def pH_to_absorbance_model_4pl(pH,upper_asymptote,Hill_slope,inflection_point,lower_asymptote):
        #This is the model that we're going to try to fit using scipy's curve_fit function.
        return lower_asymptote + (upper_asymptote - lower_asymptote) / (1 + (pH / inflection_point)**Hill_slope)

    if plot==True:
        handles_and_labels = {}

        
        def change_colors():
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
            return color_map
        color_map = change_colors()

        inflection_points_1trig = []
        inflection_points_2trig = []

        for construct_and_date in grouped_dict:
            """
            plt.fill_between(
                fitted_curve,
                fitted_curve - lower_confidence_limit,
                fitted_curve + upper_confidence_limit,
                color=color_band,
                alpha=0.45,
                label="95% CI band",
            )"""
            try:
                print(f'Construct_and_date: {construct_and_date}')
                # This is where we calculate the 95% CI at each pH for the dataset containing all of the replicates. We need to do this for each pH.
                pH_CI_band = {}
                pH_values = []
                #medians = []
                means = []
                means.append(mean)
                #median = np.median(grouped_dict[construct_and_date][pH])
                #medians.append(median)
                
                """                
                    alpha_bonferroni = 0.05/11
                    N_replicates = len(pH)
                    numerator = 0.0
                    df = 0    
                    for rep in pH:
                        n_rep = len(rep)
                        s_rep = np.std(rep, ddof=1)
                        numerator += ((n_rep - 1)*(s_rep**2))
                        df += (n_rep - 1)
                    pooled_sd = np.sqrt(numerator/df)
                    SE = pooled_sd/(np.sqrt(N_replicates))
                    t_stat = t.ppf(1-(alpha_bonferroni/2),df=df)
                    error_bar = SE*t_stat
                    pH_CI_band(pH).append(error_bar)
                    """                
                
                pH_linspace = np.linspace(min(pH_values),max(pH_values),400)
                #print(f'medians: {medians}')
                #stats.mannwhitneyu()
                #initial_guesses = [max(medians),-4,5.5,min(medians)]
                initial_guesses = [max(means),-4,5.5,min(means)]
                #The returned values from the curve_fit match the order of those passed in as initial guesses, namely: [max,slope,inflection point,min]
                best_fit_parameters,covariance_matrix = curve_fit(pH_to_absorbance_model_4pl,pH_values,means,p0=initial_guesses)
                #if "/" in construct_and_date:
                if "2Trig" in construct_and_date and best_fit_parameters[2] > 4.0:
                    inflection_points_2trig.append(best_fit_parameters[2])
                elif "2Trig" in construct_and_date and best_fit_parameters[2] < 4.0:
                    inflection_points_2trig.append(float(4.0))
                elif not "2Trig" in construct_and_date and best_fit_parameters[2] < 4.0:
                    inflection_points_1trig.append(float(4.0))
                else:
                    inflection_points_1trig.append(best_fit_parameters[2])

                #don't plot the curve_fit if the inflection point is below 4.0
                if best_fit_parameters[2] < 4.0:
                    raise ValueError(f"Calculated inflection point ({round(best_fit_parameters[2],2)}) for {construct_and_date} below 4.0")
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
                #The asterisk unpacks all of the variables in bets_fit_parameters so that the model has all the constants it needs to calculate the y value 
                #at the inflection point.
                if plot_singly:
                    if len(inflection_points_1trig)>0:
                        inflection = plt.axvline(inflection_points_1trig[-1],color="blue",linestyle="--")
                        handles_and_labels[inflection] = f"Inflection point at: {inflection_points_1trig[-1].round(2)}"
                    if len(inflection_points_2trig)>0:
                        inflection_2 = plt.axvline(inflection_points_2trig[-1],color="red",linestyle="--")
                        handles_and_labels[inflection_2] = f"Inflection point at: {inflection_points_2trig[-1].round(2)}"
                if scatter_handle is not None:
                    handles_and_labels[(line_handle, scatter_handle)] = construct_and_date
                else:
                    handles_and_labels[line_handle] = construct_and_date
            except Exception as e:
                #print(f'The curve fit for {construct_and_date} failed due to: {e}. Plotting points instead of fitted curve...')
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
                handles_and_labels[plot] = construct_and_date
            if plot_singly:
                if len(proteins_and_tests)==1:
                    plt.title(f'{protein_constructs[0]} {tests[0]}')
                    plt.ylabel(tests[0])
                    plt.xlabel('pH')
                handles = [handle for handle in handles_and_labels]
                labels = [label for handle in handles_and_labels for label in [handles_and_labels[handle]]]
                plt.legend(handles,labels,handler_map={tuple: HandlerTuple(ndivide=2, pad=1)})
                #plt.savefig(f'{construct_and_date.replace("/","_")}.png',dpi=300,bbox_inches="tight")
                plt.show()
                plt.close()
                handles_and_labels = {}

        if not plot_singly:
            if len(inflection_points_1trig)>1:
                single_trigger_inflection_point = np.mean(inflection_points_1trig)
                inflection_1 = plt.axvline(single_trigger_inflection_point,color="blue",linestyle="--")
                handles_and_labels[inflection_1] = f"Inflection point at: {single_trigger_inflection_point.round(2)}"
            elif len(inflection_points_1trig)==1:
                inflection_1 = plt.axvline(inflection_points_1trig[0],color='blue',linestyle="--")
                handles_and_labels[inflection_1] = f"Inflection point at: {np.float64(inflection_points_1trig[0]).round(2)}"

            if len(inflection_points_2trig)>1:
                double_trigger_inflection_point = np.mean(inflection_points_2trig)
                inflection_2 = plt.axvline(double_trigger_inflection_point,color="red",linestyle="--")
                handles_and_labels[inflection_2] = f"Inflection point at: {double_trigger_inflection_point.round(2)}"
            elif len(inflection_points_2trig)==1:
                inflection_2 = plt.axvline(inflection_points_2trig[0],color='red',linestyle="--")
                handles_and_labels[inflection_2] = f"Inflection point at: {np.float64(inflection_points_2trig[0]).round(2)}"        

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
            handles = [handle for handle in handles_and_labels]
            labels = [label for handle in handles_and_labels for label in [handles_and_labels[handle]]]
            plt.legend(handles,labels,handler_map={tuple: HandlerTuple(ndivide=2, pad=1)})
            #plt.savefig(f'{construct_and_date.replace("/","_")}.png',dpi=300,bbox_inches="tight")
            plt.show()
            plt.close()


#for protein_construct in protein_dict:
#    if '2Trig' not in protein_construct and 'Gravity' not in protein_construct:
        #Single-trigger A400
#        statisticize([(f'{protein_construct}','A400')],boxplot=False,plot=True,plot_singly=True)
#        statisticize([(f'2Trig-{protein_construct}','A400')],boxplot=False,plot=True,plot_singly=True)

for protein_construct in protein_dict:
#    if '2Trig' not in protein_construct and 'Gravity' not in protein_construct:
        #Single vs. Double-trigger A400
#        statisticize([(f'{protein_construct}','A400'),(f'2Trig-{protein_construct}','A400')],plot=True)
        statisticize([(f'{protein_construct}','A400'),(f'2Trig-{protein_construct}','A400')],plot=True)

#for protein_construct in protein_dict:
 #   if '2Trig' not in protein_construct and 'Gravity' not in protein_construct:
        #Single vs. Double-trigger A280-1hr
  #      statisticize([(f'{protein_construct}','A280_1hr'),(f'2Trig-{protein_construct}','A280_1hr')],boxplot=False,plot=True)
    #    statisticize([(f'{protein_construct}','A280_1hr'),(f'2Trig-{protein_construct}','A280_1hr')],boxplot=False,plot=True)

#for protein_construct in protein_dict:
    #if '2Trig' not in protein_construct and 'Gravity' not in protein_construct:
        #Single vs. Double-trigger A280_48-72hr
        #statisticize([(f'{protein_construct}','A280_48-72hr'),(f'2Trig-{protein_construct}','A280_48-72hr')],boxplot=False,plot=True)
        #statisticize([(f'{protein_construct}','A280_48-72hr'),(f'2Trig-{protein_construct}','A280_48-72hr')],boxplot=False,plot=True)

#for protein_construct in protein_dict:
 #   if '2Trig' not in protein_construct and 'Gravity' not in protein_construct:
  #      #Single trigger A280_1hr vs A280_48-72hr
   #     statisticize([(f'{protein_construct}','A280_1hr'),(f'{protein_construct}','A280_48-72hr')],boxplot=False,plot=True,require_matched_groups_var=True)

#Single trigger TV-vWA A400
#statisticize([('10xHis-1TEL-TV-vWA','A400'),('10xHis-1TEL-TV-vWA (Gravity)','A400')],plot=True)
#statisticize([('10xHis-1TEL-TV-vWA','A400'),('10xHis-1TEL-TV-vWA (Gravity)','A400')],plot=True)
