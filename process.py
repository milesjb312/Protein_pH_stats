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

def format_protein_title(construct_name: str) -> str:
    name = construct_name.strip()

    # Remove prefixes that should not appear in the title
    name = name.replace("2Trig-", "")
    name = name.replace(" (Gravity)", "")

    parts = name.split("-")

    # 1) tag
    tag = "10xHis-" if parts[0] == "10xHis" else ""

    # 2) scaffold
    # You said this will always be 1TEL
    scaffold = "1TEL"

    # 3) linker
    # Expected examples: GG, PA, SR, TV
    linker = ""
    if "GG" in parts:
        linker = "GG"
    elif "PA" in parts:
        linker = "PA"
    elif "SR" in parts:
        linker = "[sr]"
    elif "TV" in parts:
        linker = "TV"

    # 4) target protein
    if "TNK1.UBA" in name:
        target_protein = "TNK1.UBA"
    elif "DARPin" in name:
        target_protein = "DARPin"
    elif "vWA" in name:
        target_protein = "vWA"
    else:
        target_protein = "Unknown"

    return f"{tag}{scaffold}-{linker}-{target_protein}"

def format_legend_label(construct, test_name, protein_constructs, tests, n):

    is_2trig = "2Trig" in construct
    is_gravity = "Gravity" in construct
    is_vwa = any(("vWA" in c) and ("2Trig" not in c) for c in protein_constructs)

    has_1hr = "A280_1hr" in tests
    has_48hr = "A280_48-72hr" in tests

    # Case 1: 1Trig vs 2Trig
    if len(protein_constructs) == 2 and len(tests) == 1:
        if is_vwa and any("Gravity" in c for c in protein_constructs):
            base_label = "Gravity" if is_gravity else "ÄKTA"

            if test_name == "A280_1hr":
                label = f"{base_label} 1 Hour"
            elif test_name == "A280_48-72hr":
                label = f"{base_label} 48 Hour"
            else:
                label = base_label

        else:
            label = "2Trig" if is_2trig else "1Trig"
    # Case 2: 1 hr vs 48 hr
    elif len(protein_constructs) == 1 and has_1hr and has_48hr:
        trig_label = "2Trig" if is_2trig else "1Trig"
        time_label = "1 Hour" if test_name == "A280_1hr" else "48 Hour"
        label = f"{trig_label} {time_label}"
    else:
        label = construct
    return f"{label} (n = {n})"

def statisticize(proteins_and_tests:list,plot_singly=False,plot=False,tailored=False,pool=False):
    """
    Normally, you should pass only one test and one protein construct into 'proteins_and_tests', in the form of a list of tuples, where the
    tuple looks like: (protein_construct,test)
    'protein_constructs' are the protein_dict dictionary keys that refer to the particular protein constructs of interest.
    'tests' are the types of tests you performed, referred to by the variable being measured, or the column name in the csv. Ie.: 'A400', 'A280_1hr'    
    """
    #The first half of this function is just a filter pulling applicable data from the protein_dict dictionary.
    data = {'protein construct':[],'date':[],'pH':[],'test':[],'data':[],'tailored_data':[], 'is_starred':[]}
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
                            value = (theoretical_max_concentration-replicate)/theoretical_max_concentration
                        else:
                            value = replicate/theoretical_max_concentration
                        data['tailored_data'].append(value)
                        data['data'].append(value)
                        data['is_starred'].append(False)
                    else:
                        if test in ['A280_1hr','A280_48-72hr']:
                            #For precipitation assays, find the amount of protein left as a percentage.
                            value = (theoretical_max_concentration-float(replicate.strip("*")))/theoretical_max_concentration
                        else:
                            value = float(replicate.strip("*"))/theoretical_max_concentration
                        data['data'].append(value)
                        data['tailored_data'].append(value)   # keep numeric
                        data['is_starred'].append(True)

    #--------------------Now that the main filtering is done, it's time to make confidence intervals and fit curves for each set of construct:pH values---------------#

    pHs = []
    for pH in data['pH']:
        if pH not in pHs:
            pHs.append(pH)

    if tailored == True:
        data_quality = 'tailored_data'
    else:
        data_quality = 'data'

    constructs_by_pH_by_date = {}
    for replicate in range(len(data['protein construct'])):
        construct = data['protein construct'][replicate]
        test_name = data['test'][replicate]
        construct_key = (construct, test_name)
        pH = data['pH'][replicate]
        date = data['date'][replicate]
        measurement = data[data_quality][replicate]
        if construct_key not in constructs_by_pH_by_date:
            constructs_by_pH_by_date[construct_key] = {}
        is_starred = data['is_starred'][replicate]
        if pH not in constructs_by_pH_by_date[construct_key]:
            constructs_by_pH_by_date[construct_key][pH] = {"dates": {}, "dates_starred": {}}
        if date not in constructs_by_pH_by_date[construct_key][pH]["dates"]:
            constructs_by_pH_by_date[construct_key][pH]["dates"][date] = []
            constructs_by_pH_by_date[construct_key][pH]["dates_starred"][date] = []
        constructs_by_pH_by_date[construct_key][pH]["dates"][date].append(measurement)
        constructs_by_pH_by_date[construct_key][pH]["dates_starred"][date].append(is_starred)

    biological_n_by_construct = {}

    for construct in constructs_by_pH_by_date:
        unique_dates = set()
        for pH in constructs_by_pH_by_date[construct]:
            for date in constructs_by_pH_by_date[construct][pH]["dates"]:
                unique_dates.add(date)
        biological_n_by_construct[construct] = len(unique_dates)
      
    # Compute the 95% CI Band for each pH value within each construct
    alpha_bonferroni = 0.05/11
    for construct in constructs_by_pH_by_date:
        for pH in constructs_by_pH_by_date[construct]:
            # Eli Note: The line below gives back a list of lists. Each outer list is a biological replicate and each inner list includes the measurements from the nanodrop for that biological replicate.
            pH_dates = constructs_by_pH_by_date[construct][pH]["dates"]
            pH_starred = constructs_by_pH_by_date[construct][pH]["dates_starred"]
            if pool:
                weighted_variances = 0.0
                pooled_df = 0
            all_drops_at_current_pH = []
            for date_key in pH_dates:
                drops = pH_dates[date_key]
                star_flags = pH_starred[date_key]
                if tailored:
                    drops = [d for d, s in zip(drops, star_flags) if not s]
                if len(drops) == 0:
                    continue
                drops = np.array(drops, dtype=float)
                all_drops_at_current_pH.extend(drops)
                if pool:
                    sd_of_drops = np.std(drops, ddof=1)
                    weighted_variances += ((len(drops) - 1)*(sd_of_drops**2))
                    pooled_df += (len(drops) - 1)
            if len(all_drops_at_current_pH) == 0:
                continue
            N_replicates = len(all_drops_at_current_pH)
            if pool:
                pooled_sd = np.sqrt(weighted_variances/pooled_df)
                SE = pooled_sd/(np.sqrt(N_replicates))
                t_stat = t.ppf(1-(alpha_bonferroni/2),df=pooled_df)
            else:
                sd = np.std(all_drops_at_current_pH)
                SE = sd/(np.sqrt(N_replicates))
                t_stat = t.ppf(1-(alpha_bonferroni/2),df=len(all_drops_at_current_pH)-1)
            print(f'SE: {SE}')
            error_bar = SE*t_stat
            mean_value = np.mean(all_drops_at_current_pH)
            constructs_by_pH_by_date[construct][pH]["mean"] = mean_value
            if pool:
                constructs_by_pH_by_date[construct][pH]["pooled_sd"] = pooled_sd
                constructs_by_pH_by_date[construct][pH]["df"] = pooled_df
            else:
                constructs_by_pH_by_date[construct][pH]["sd"] = sd
                constructs_by_pH_by_date[construct][pH]['df'] = len(all_drops_at_current_pH)-1
            constructs_by_pH_by_date[construct][pH]["SE"] = SE
            constructs_by_pH_by_date[construct][pH]["t_stat"] = t_stat
            constructs_by_pH_by_date[construct][pH]["error_bar"] = error_bar
            constructs_by_pH_by_date[construct][pH]["CI_lower"] = mean_value - error_bar
            constructs_by_pH_by_date[construct][pH]["CI_upper"] = mean_value + error_bar
            constructs_by_pH_by_date[construct][pH]["N"] = N_replicates
        
    def pH_to_absorbance_model_4pl(pH,upper_asymptote,Hill_slope,inflection_point,lower_asymptote):
        #This is the model (4 Parameter Logistic, 4PL) that we're going to try to fit using scipy's curve_fit function.
        return lower_asymptote + (upper_asymptote - lower_asymptote) / (1 + (pH / inflection_point)**Hill_slope)

    if plot==True:
        # Make legend smaller:
            # GG A400 1 vs 2
            # SR 1 HR 1 vs 2
            # GG 1 HR 1 vs 2
            # vWA 1 HR non-gravity vs gravity
        # Legend in bottom:
            # Why 2Trig-PA in bottom by not 1Trig 1HR vs 48 HR?

        plt.figure(figsize=(8.75,6))
        handles_and_labels = {}

        def change_colors():
            keys_1trig = []
            keys_2trig = []
            keys_48_hour = []
            keys_2Trig_48_hour = []
            keys_Gravity = []
            for construct_key in constructs_by_pH_by_date:
                construct, test_name = construct_key
                if "Gravity" in construct:
                    keys_Gravity.append(construct_key)
                elif test_name == "A280_48-72hr" and not '2Trig' in construct:
                    keys_48_hour.append(construct_key)
                elif test_name == 'A280_48-72hr' and '2Trig' in construct:
                    keys_2Trig_48_hour.append(construct_key)
                elif "2Trig" in construct:
                    keys_2trig.append(construct_key)
                else:
                    keys_1trig.append(construct_key)
            def shade_list(cmap, n, lo=0.25, hi=0.9):
                if n<= 1:
                    return [cmap((lo + hi) / 2)]
                return [cmap(x) for x in np.linspace(lo, hi, n)]
            blue_shades = shade_list(plt.cm.Blues, len(keys_1trig))
            orange_shades = shade_list(plt.cm.Oranges, len(keys_2trig))
            crimson_shades = shade_list(plt.cm.PuRd, len(keys_48_hour))
            purple_shades = shade_list(plt.cm.Purples, len(keys_2Trig_48_hour))
            lime_shades = shade_list(plt.cm.YlGn, len(keys_Gravity))
            color_map = {}
            for construct, color in zip(sorted(keys_1trig), blue_shades):
                color_map[construct] = color
            for construct, color in zip(sorted(keys_2trig), orange_shades):
                color_map[construct] = color
            for construct, color in zip(sorted(keys_48_hour), crimson_shades):
                color_map[construct] = color
            for construct, color in zip(sorted(keys_2Trig_48_hour), purple_shades):
                color_map[construct] = color
            for construct, color in zip(sorted(keys_Gravity), lime_shades):
                color_map[construct] = color
            return color_map
        color_map = change_colors()

        inflection_points_1trig_1hr = []
        inflection_points_1trig_48hr = []
        inflection_points_2trig_1hr = []
        inflection_points_2trig_48hr = []
        inflection_points_gravity = []
        any_starred_points = False

        for construct_key in constructs_by_pH_by_date:
            construct, test_name = construct_key
            pH_values = []
            means = []
            ci_lower_bounds = []
            ci_upper_bounds = []
            construct_pHs =[]
            construct_data = []
            construct_starred = []
            # sort pHs numerically
            sorted_pHs = sorted(constructs_by_pH_by_date[construct_key].keys(), key=float)
            for pH in sorted_pHs:
                pH_entry = constructs_by_pH_by_date[construct_key][pH]
                # Always keep raw points for plotting
                for date, measurements in pH_entry["dates"].items():
                    star_flags = pH_entry["dates_starred"][date]
                    for measurement, is_starred in zip(measurements, star_flags):
                        construct_pHs.append(float(pH))
                        construct_data.append(measurement)
                        construct_starred.append(is_starred)
                # Only use pHs with computed statistics for fitting / CI bands
                if "mean" not in pH_entry:
                    continue
                pH_values.append(float(pH))
                means.append(pH_entry["mean"])
                ci_lower_bounds.append(pH_entry["CI_lower"])
                ci_upper_bounds.append(pH_entry["CI_upper"])

            try:
                pH_linspace = np.linspace(min(pH_values),max(pH_values),400)
                initial_guesses = [max(means),-4,5.5,min(means)]
                #The returned values from the curve_fit match the order of those passed in as initial guesses, namely: [max,slope,inflection point,min]
                best_fit_parameters,covariance_matrix = curve_fit(pH_to_absorbance_model_4pl,pH_values,means,p0=initial_guesses)
                print(best_fit_parameters)
                inflection_value = best_fit_parameters[2] if best_fit_parameters[2] >= 4.0 else 4.0
                if "Gravity" in construct:
                    inflection_points_gravity.append(inflection_value)
                elif "2Trig" in construct:
                    if test_name == "A280_48-72hr":
                        inflection_points_2trig_48hr.append(inflection_value)
                    else:
                        inflection_points_2trig_1hr.append(inflection_value)
                else:
                    if test_name == "A280_48-72hr":
                        inflection_points_1trig_48hr.append(inflection_value)
                    else:
                        inflection_points_1trig_1hr.append(inflection_value)

                ###-------------------------------###
                #To plot the fitted curve:
                #don't plot the curve_fit if the inflection point is below 4.0
                if best_fit_parameters[2] < 4.0:
                    raise ValueError(f"Calculated inflection point ({round(best_fit_parameters[2],2)}) for {construct} below 4.0")
                
                plot = plt.plot(pH_linspace,pH_to_absorbance_model_4pl(pH_linspace,*best_fit_parameters), color=color_map[construct_key])

                line_handle = plot[0]

                ###------plot CI band around the mean values at each pH---------------------###
                if "Gravity" in construct:
                    plt.fill_between(pH_values, ci_lower_bounds, ci_upper_bounds, color='lime', alpha=0.20)
                elif test_name == "A280_48-72hr" and "2Trig" not in construct:
                    plt.fill_between(pH_values, ci_lower_bounds, ci_upper_bounds, color='pink', alpha=0.20)
                elif test_name == "A280_48-72hr" and "2Trig" in construct:
                    plt.fill_between(pH_values, ci_lower_bounds, ci_upper_bounds, color='mediumpurple', alpha=0.20)
                elif "2Trig" in construct:
                    plt.fill_between(pH_values, ci_lower_bounds, ci_upper_bounds, color='tomato', alpha=0.20)
                else:
                    plt.fill_between(pH_values, ci_lower_bounds, ci_upper_bounds, color='cyan', alpha=0.20)

                # plot raw scatter points
                jitter_strength = 0.08
                jittered_pHs = np.array(construct_pHs) + np.random.uniform(-jitter_strength, jitter_strength, len(construct_pHs))
                if "2Trig" in construct and test_name == 'A280_48-72hr':
                    scatter_handle = plt.scatter(jittered_pHs,construct_data,marker="^", color=color_map[construct_key], alpha=0.3)
                elif "2Trig" in construct:
                    scatter_handle = plt.scatter(jittered_pHs,construct_data,marker="^", color=color_map[construct_key], alpha=0.2)
                else:
                    scatter_handle = plt.scatter(jittered_pHs,construct_data, color=color_map[construct_key], alpha=0.2)
                # Plot asterisks for omitted (starred) data
                starred_jittered_pHs = [x for x, s in zip(jittered_pHs, construct_starred) if s]
                starred_data = [y for y, s in zip(construct_data, construct_starred) if s]
                if len(starred_data) > 0:
                    for x, y in zip(starred_jittered_pHs, starred_data):
                        plt.annotate("*",(x, y),xytext=(0.005, 0.00000005),textcoords="offset points",ha="left",va="bottom",color="black",fontsize=11,fontweight="light",zorder=6)
                    any_starred_points = True

                # mean points with Bonferroni error bars
                #plt.errorbar(x=pH_values,y=means, yerr=[np.array(means) - np.array(ci_lower_bounds), np.array(ci_upper_bounds) - np.array(means)], fmt='none', ecolor=color_map[construct], capsize=4, alpha=0.9)
                ###--------------------------------------------------------------------------###

                if plot_singly:
                    if len(inflection_points_1trig)>0:
                        inflection = plt.axvline(inflection_points_1trig[-1],color="blue",linestyle="--")
                        handles_and_labels[inflection] = f"Inflection point: pH = {inflection_points_1trig[-1].round(2)}"
                    if len(inflection_points_2trig)>0:
                        inflection_2 = plt.axvline(inflection_points_2trig[-1],color="red",linestyle="--")
                        handles_and_labels[inflection_2] = f"Inflection point: pH = {inflection_points_2trig[-1].round(2)}"
                legend_label = format_legend_label(construct,test_name,protein_constructs,tests,biological_n_by_construct[construct_key])
                handles_and_labels[(line_handle, scatter_handle)] = legend_label

            except Exception as e:
                print(f'The curve fit for {construct} failed due to: {e}. Plotting points instead of fitted curve...')
                construct_pHs = []
                construct_data = []
                construct_starred = []
                sorted_pHs = sorted(constructs_by_pH_by_date[construct_key].keys(), key=float)
                for pH in sorted_pHs:
                    pH_entry = constructs_by_pH_by_date[construct_key][pH]
                    for date, measurements in pH_entry["dates"].items():
                        star_flags = pH_entry["dates_starred"][date]
                        for measurement, is_starred in zip(measurements, star_flags):
                            construct_pHs.append(float(pH))
                            construct_data.append(measurement)
                            construct_starred.append(is_starred)
                if "Gravity" in construct:
                    plt.fill_between(pH_values, ci_lower_bounds, ci_upper_bounds, color='lime', alpha=0.20)
                elif test_name == "A280_48-72hr" and "2Trig" not in construct:
                    plt.fill_between(pH_values, ci_lower_bounds, ci_upper_bounds, color='pink', alpha=0.20)
                elif test_name == "A280_48-72hr" and "2Trig" in construct:
                    plt.fill_between(pH_values, ci_lower_bounds, ci_upper_bounds, color='mediumpurple', alpha=0.20)
                elif "2Trig" in construct:
                    plt.fill_between(pH_values, ci_lower_bounds, ci_upper_bounds, color='tomato', alpha=0.20)
                else:
                    plt.fill_between(pH_values, ci_lower_bounds, ci_upper_bounds, color='cyan', alpha=0.20)
                jitter_strength = 0.08
                jittered_pHs = np.array(construct_pHs) + np.random.uniform(-jitter_strength, jitter_strength, len(construct_pHs))
                if "2Trig" in construct and test_name == 'A280_48-72hr':
                    scatter_handle = plt.scatter(jittered_pHs,construct_data,marker="^", color=color_map[construct_key], alpha=0.3)
                elif "2Trig" in construct:
                    scatter_handle = plt.scatter(jittered_pHs,construct_data,marker="^", color=color_map[construct_key], alpha=0.2)
                else:
                    scatter_handle = plt.scatter(jittered_pHs,construct_data, color=color_map[construct_key], alpha=0.2)
                starred_jittered_pHs = [x for x, s in zip(jittered_pHs, construct_starred) if s]
                starred_data = [y for y, s in zip(construct_data, construct_starred) if s]
                if len(starred_data) > 0:
                    for x, y in zip(starred_jittered_pHs, starred_data):
                        plt.annotate("*",(x, y),xytext=(0.005, 0.00000005),textcoords="offset points",ha="left",va="bottom",color="black",fontsize=11,fontweight="light",zorder=6)
                    any_starred_points = True
                legend_label = format_legend_label(construct,test_name,protein_constructs,tests,biological_n_by_construct[construct_key])
                handles_and_labels[scatter_handle] = legend_label

            # mean points with Bonferroni error bars
            #plt.errorbar(x=pH_values,y=means, yerr=[np.array(means) - np.array(ci_lower_bounds), np.array(ci_upper_bounds) - np.array(means)], fmt='none', ecolor=color_map[construct], capsize=4, alpha=0.9)
            ###--------------------------------------------------------------------------###

            if plot_singly:
                if len(proteins_and_tests)==1:
                    plt.title(f'{protein_constructs[0]} {tests[0]}')
                    plt.ylabel(tests[0])
                    plt.xlabel('pH')
                handles = [handle for handle in handles_and_labels]
                labels = [handles_and_labels[handle] for handle in handles_and_labels]
                plt.legend(handles,labels,handler_map={tuple: HandlerTuple(ndivide=2, pad=1)},loc="center left",bbox_to_anchor=(1.02,0.5))
                #plt.savefig(f'{construct_and_date.replace("/","_")}.png',dpi=300,bbox_inches="tight")
                plt.show()
                plt.close()
                handles_and_labels = {}

        if not plot_singly:
            is_gravity_comparison = (len(protein_constructs) == 2 and len(tests) == 1 and any("Gravity" in construct for construct in protein_constructs))
            is_time_comparison = (len(protein_constructs) == 1 and "A280_1hr" in tests and "A280_48-72hr" in tests)
            if is_time_comparison:
                if "2Trig" in protein_constructs[0]:
                    if len(inflection_points_2trig_1hr) >= 1:
                        inflection_1hr = plt.axvline(np.mean(inflection_points_2trig_1hr), color="darkorange", linestyle="--")
                        handles_and_labels[inflection_1hr] = f"1 Hour Inflection: pH = {np.mean(inflection_points_2trig_1hr):.2f}"
                    if len(inflection_points_2trig_48hr) >= 1:
                        inflection_48hr = plt.axvline(np.mean(inflection_points_2trig_48hr), color="indigo", linestyle="--")
                        handles_and_labels[inflection_48hr] = f"48 Hour Inflection: pH = {np.mean(inflection_points_2trig_48hr):.2f}"
                else:
                    if len(inflection_points_1trig_1hr) >= 1:
                        inflection_1hr = plt.axvline(np.mean(inflection_points_1trig_1hr), color="blue", linestyle="--")
                        handles_and_labels[inflection_1hr] = f"1 Hour Inflection: pH = {np.mean(inflection_points_1trig_1hr):.2f}"
                    if len(inflection_points_1trig_48hr) >= 1:
                        inflection_48hr = plt.axvline(np.mean(inflection_points_1trig_48hr), color="purple", linestyle="--")
                        handles_and_labels[inflection_48hr] = f"48 Hour Inflection: pH = {np.mean(inflection_points_1trig_48hr):.2f}"
            elif is_gravity_comparison:
                # non-gravity line
                if "A400" in tests:
                    nongravity_points = inflection_points_1trig_1hr
                    nongravity_color = "blue"
                    nongravity_label = "ÄKTA Inflection"
                elif "A280_48-72hr" in tests:
                    nongravity_points = inflection_points_1trig_48hr
                    nongravity_color = "purple"
                    nongravity_label = "ÄKTA 48 Hour Inflection"
                elif "A280_1hr" in tests:
                    nongravity_points = inflection_points_1trig_1hr
                    nongravity_color = "blue"
                    nongravity_label = "ÄKTA 1 Hour Inflection"
                if len(nongravity_points) >= 1:
                    inflection_ng = plt.axvline(np.mean(nongravity_points), color=nongravity_color, linestyle="--")
                    handles_and_labels[inflection_ng] = f"{nongravity_label}: pH = {np.mean(nongravity_points):.2f}"

                if len(inflection_points_gravity) >= 1:
                    inflection_g = plt.axvline(np.mean(inflection_points_gravity), color="green", linestyle="--")
                    if "A400" in tests:
                        handles_and_labels[inflection_g] = f"Gravity Inflection: pH = {np.mean(inflection_points_gravity):.2f}"
                    elif "A280_48-72hr" in tests:
                        handles_and_labels[inflection_g] = f"Gravity 48 Hour Inflection: pH = {np.mean(inflection_points_gravity):.2f}"
                    elif "A280_1hr" in tests:
                        handles_and_labels[inflection_g] = f"Gravity 1 Hour Inflection: pH = {np.mean(inflection_points_gravity):.2f}"
            else:
                # standard 1Trig vs 2Trig comparison
                if len(inflection_points_1trig_1hr) >= 1:
                    if "A280_48-72hr" in tests:
                        color_1trig = "purple"
                    else:
                        color_1trig = "blue"
                    inflection_1 = plt.axvline(np.mean(inflection_points_1trig_1hr), color=color_1trig, linestyle="--")
                    handles_and_labels[inflection_1] = f"1Trig Inflection: pH = {np.mean(inflection_points_1trig_1hr):.2f}"
                if len(inflection_points_2trig_1hr) >= 1:
                    if "A280_48-72hr" in tests:
                        color_2trig = "indigo"
                    else:
                        color_2trig = "red"
                    inflection_2 = plt.axvline(np.mean(inflection_points_2trig_1hr), color=color_2trig, linestyle="--")
                    handles_and_labels[inflection_2] = f"2Trig Inflection: pH = {np.mean(inflection_points_2trig_1hr):.2f}" 
            if "A400" == tests[0]:
                y_axis = 'Absorbance'
                assay = 'Turbidity'
            else:
                y_axis = "Precipitation Loss / mg/mL"
                assay = 'Solubility'
            if len(protein_constructs) == 1:
                protein_title = format_protein_title(protein_constructs[0])
                plt.title(f"{protein_title}, {assay}", fontsize = 14)
                plt.ylabel(y_axis)
                plt.xlabel("pH")
            else:
                if len(protein_constructs) == 2 and len(tests) == 1:
                    # Use the first construct as the base protein title
                    protein_title = format_protein_title(protein_constructs[0])
                    plt.title(f"{protein_title}, {assay}", fontsize = 14)
                    plt.ylabel(y_axis)
                    plt.xlabel("pH")
                elif len(protein_constructs) == 1 and len(tests) == 2:
                    protein_title = format_protein_title(protein_constructs[0])
                    plt.title(f"{protein_title}, {assay}", fontsize = 14)
                    plt.ylabel(y_axis)
                    plt.xlabel("pH")
            # add asterisk legend entry LAST so it appears on bottom row
            if any_starred_points:
                    asterisk_handle = plt.scatter([], [], marker='$*$', color='black', s=40)
                    handles_and_labels[asterisk_handle] = "Omitted from curve fitting"
            handles = [handle for handle in handles_and_labels]
            labels = [label for handle in handles_and_labels for label in [handles_and_labels[handle]]]
            is_pa = any("PA-TNK1.UBA" in construct for construct in protein_constructs)
            is_vwa = any("vWA" in construct for construct in protein_constructs)
            is_1hr_only = ("A280_1hr" in tests and "A280_48-72hr" not in tests)
            is_1hr_vs_48hr = ("A280_1hr" in tests and "A280_48-72hr" in tests)
            is_pa_1trig_vs_2trig_1hr = is_pa and len(protein_constructs) == 2 and is_1hr_only
            is_vwa_1trig_vs_2trig_1hr = is_vwa and len(protein_constructs) == 2 and is_1hr_only
            is_pa_1trig_1hr_vs_48hr = is_pa and len(protein_constructs) == 1 and "2Trig" not in protein_constructs[0] and is_1hr_vs_48hr
            is_vwa_1trig_1hr_vs_48hr = (is_vwa and len(protein_constructs) == 1 and "2Trig" not in protein_constructs[0] and is_1hr_vs_48hr)
            is_vwa_2trig_1hr_vs_48hr = (is_vwa and len(protein_constructs) == 1 and "2Trig" in protein_constructs[0] and is_1hr_vs_48hr)
            is_vwa_1hr_non_gravity_vs_gravity = (is_vwa and len(protein_constructs) == 2 and is_1hr_only and any("Gravity" in construct for construct in protein_constructs))
            
            if (is_pa_1trig_vs_2trig_1hr or is_vwa_1trig_vs_2trig_1hr or is_pa_1trig_1hr_vs_48hr or is_vwa_1trig_1hr_vs_48hr or is_vwa_2trig_1hr_vs_48hr or is_vwa_1hr_non_gravity_vs_gravity):
                plt.legend(handles, labels,handler_map={tuple: HandlerTuple(ndivide=2, pad=1)},loc="lower right")
            else:
                plt.legend(handles, labels,handler_map={tuple: HandlerTuple(ndivide=2, pad=1)},loc="upper right")
            plt.savefig(f'{protein_construct.replace("/","_")}_{test}.png',dpi=300,bbox_inches="tight")
            plt.show()
            plt.close()


for protein_construct in protein_dict:
    if '2Trig' not in protein_construct and 'Gravity' not in protein_construct:
        #Single-trigger A400
        #statisticize([(f'{protein_construct}','A400')],plot=True,plot_singly=True)
        #Double-trigger A400
        #statisticize([(f'2Trig-{protein_construct}','A400')],plot=True,plot_singly=True)
        #Single vs. Double-trigger A400
        statisticize([(f'{protein_construct}','A400'),(f'2Trig-{protein_construct}','A400')],plot=True)

for protein_construct in protein_dict:
    if '2Trig' not in protein_construct and 'Gravity' not in protein_construct:
        #Single-trigger A280
        #statisticize([(f'{protein_construct}','A280_1hr')],plot=True,plot_singly=True)
        #Double-trigger A280
        #statisticize([(f'2Trig-{protein_construct}','A280_1hr')],plot=True,plot_singly=True)
        #Single vs. Double-trigger A280-1hr
        statisticize([(f'{protein_construct}','A280_1hr'),(f'2Trig-{protein_construct}','A280_1hr')],plot=True)

for protein_construct in protein_dict:
    if '2Trig' not in protein_construct and 'Gravity' not in protein_construct:
        #Single-trigger A280_48-72hr
        #statisticize([(f'{protein_construct}','A280_48-72hr')],plot=True,plot_singly=True)
        #Double-trigger A280_48-72hr
        #statisticize([(f'2Trig-{protein_construct}','A280_48-72hr')],plot=True,plot_singly=True)
        #Single vs. Double-trigger A280_48-72hr
        #statisticize([(f'{protein_construct}','A280_48-72hr'),(f'2Trig-{protein_construct}','A280_48-72hr')],plot=True)
        #Single-trigger A280 1hr vs 48-72 hr
        statisticize([(f'{protein_construct}','A280_1hr'),(f'{protein_construct}','A280_48-72hr')],plot=True,tailored=True)
    elif '2Trig' in protein_construct:
        #Double-trigger A280 1hr vs 48-72 hr
        statisticize([(f'{protein_construct}','A280_1hr'),(f'{protein_construct}','A280_48-72hr')],plot=True,tailored=True)

#Single trigger TV-vWA A400, A280 1 HR, A280 48-72 HR
statisticize([('10xHis-1TEL-TV-vWA','A400'),('10xHis-1TEL-TV-vWA (Gravity)','A400')],plot=True)
statisticize([('10xHis-1TEL-TV-vWA','A280_1hr'),('10xHis-1TEL-TV-vWA (Gravity)','A280_1hr')],plot=True)
statisticize([('10xHis-1TEL-TV-vWA','A280_48-72hr'),('10xHis-1TEL-TV-vWA (Gravity)','A280_48-72hr')],plot=True)