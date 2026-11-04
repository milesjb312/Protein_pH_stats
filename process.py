# Disclosure Statement: The authors of this code used ChatGPT services to develop and edit most of the code. The authors reviewed every line of
# AI-produced code before implementation to ensure that it performed the desired functions without inappropriately changing or removing data.

import csv
from scipy.optimize import curve_fit
from scipy.stats import t
import matplotlib.pyplot as plt
from matplotlib.legend_handler import HandlerTuple
import matplotlib.ticker as mticker
import numpy as np
from pathlib import Path

protein_dict = {}

BASE_DIR = Path(__file__).resolve().parent
CSV_PATH = BASE_DIR / "data.csv"

"""
Populate our protein_dict dictionary with the specified constructs and their corresponding specified data from the csv by calling the master function 'statisticize'.
For example, we may use the command statisticize([('Construct1','A400'),('Construct2','A400')]) to populate protein_dict with turbidity data corresponding to both
single- and double-trigger constructs. Each protein construct within protein_dict will include information about each biological replicate's data, stock protein concentration
used to perform the assay, the assay type, the pH the measurement was taken at, the raw measurement, the stock ID (if applicable), and the vWA ID (if applicable). We do 
include measurements omitted from curve fitting (marked in the csv file with a *, ^, or ~) within protein_dict. To clarify, this code uses a variety of labels for single- and
double-trigger constructs, but annotations will refer to single- and double-trigger as 1Trig and 2Trig, respectively. 
"""
with open(CSV_PATH, newline="", encoding="utf-8-sig") as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        protein_construct = row['Protein Construct'].strip()
        if protein_construct not in protein_dict:
            protein_dict[protein_construct] = {}
        bio_rep_key = row['Date'] + "_" + row['Stock Concentration mg/mL'] + '_mg/mL_pH_' + row['pH']
        if bio_rep_key not in protein_dict[protein_construct]:
            protein_dict[protein_construct][bio_rep_key] = {'A400': [],'A280_1hr': [],'A280_48hr': [],'A280_72hr': [],'Stock ID': row.get('Stock ID', '').strip(),'Behavior ID': row.get('Behavior ID', '').strip()}
        if row['A400'] != "":
            raw_value = row['A400'].strip()
            if '*' in raw_value or '^' in raw_value or '~' in raw_value:
                protein_dict[protein_construct][row['Date']+"_"+row['Stock Concentration mg/mL']+'_mg/mL_pH_'+row['pH']]['A400'].append(raw_value)
            else:
                protein_dict[protein_construct][row['Date']+"_"+row['Stock Concentration mg/mL']+'_mg/mL_pH_'+row['pH']]['A400'].append(float(raw_value))
        if row['A280_1hr']!="":
            raw_value = row['A280_1hr'].strip()
            if '*' in raw_value or '^' in raw_value or '~' in raw_value:
                protein_dict[protein_construct][row['Date']+"_"+row['Stock Concentration mg/mL']+'_mg/mL_pH_'+row['pH']]['A280_1hr'].append(raw_value)
            else:
                protein_dict[protein_construct][row['Date']+"_"+row['Stock Concentration mg/mL']+'_mg/mL_pH_'+row['pH']]['A280_1hr'].append(float(raw_value))
        if row['A280_48-72hr'] != "":
            raw_value = row['A280_48-72hr'].strip()
            timepoint = row.get('A280_48-72hr Time', '').strip().lower()
            if timepoint == "72":
                time_test = "A280_72hr"
            else:
                time_test = "A280_48hr"
            if '*' in raw_value or '^' in raw_value or '~' in raw_value:
                protein_dict[protein_construct][bio_rep_key][time_test].append(raw_value)
            else:
                protein_dict[protein_construct][bio_rep_key][time_test].append(float(raw_value))

"""
The next two helper functions format the title and legend for each plot using the information in protein_dict. Each plot is titled with the base protein and the legend specifies which
points are 1Trig or 2Trig, if they are from a 1 hour or 48 hour solubility assay, and if the stock protein concentration was made using an ÄKTA or gravity nickel affinity
purification coloumn. For example, for the turbidity graph comparing the 1Trig and 2Trig variants of 1TEL-GG-TNK1.UBA, the title will say "1TEL-GG-TNK1.UBA" and the 
legend will specifiy which dots, inflection point line, 95% confidence intervals, and four parameter curves correlate to 1Trig and 2Trig.
"""
def format_chart_title(construct_name: str,multiplicity: int) -> str:
    name = construct_name.strip()
    # Remove prefixes that should not appear in the title
    #name = name.replace("2Trig-", "")
    name = name.replace(" (Gravity)", "")
    parts = name.split("-")
    # 1) tag
    #tag = "10xHis-" if parts[0] == "10xHis" else ""
    # 2) scaffold
    # You said this will always be 1TEL
    scaffold = "1TEL"
    # 3) linker
    # Expected examples: GG, PA, SR, TV
    linker = ""
    if "2Trig" in parts:
        trigger = "dt"
    else:
        trigger = "st"
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
    elif "Knockout" in name:
        target_protein = 'Knockout'
    elif 'Alone' in name:
        target_protein = 'Alone'
    else:
        target_protein = "Unknown"
    if multiplicity == 1:
        return f"{trigger}-{scaffold}-{linker}-{target_protein}"
    else:
        return f"{scaffold}-{linker}-{target_protein}"

def time_label_from_test(test_name):
    if test_name == "A280_1hr":
        return "1 Hour"
    elif test_name == "A280_48hr":
        return "48 Hour"
    elif test_name == "A280_72hr":
        return "72 Hour"
    return ""

def format_legend_label(construct, test_name, protein_constructs, tests, n):
    is_2trig = "2Trig" in construct
    is_gravity = "Gravity" in construct
    is_vwa_1trig = any(("vWA" in c) and ("2Trig" not in c) for c in protein_constructs)
    is_alone_knockout_comparison = {
        "1TEL Alone", "2Trig-1TEL Alone", "1TEL Knockout"
    }.issubset(protein_constructs)
    has_1hr = "A280_1hr" in tests
    has_later_timepoint = any(t in tests for t in ["A280_48hr", "A280_72hr"])
    is_solubility = test_name in ["A280_1hr", "A280_48hr", "A280_72hr"]
    # Case 1: comparing two constructs with one assay/timepoint
    if (len(protein_constructs) == 2 or is_alone_knockout_comparison) and len(tests) == 1:
        if is_vwa_1trig and any("Gravity" in c for c in protein_constructs):
            base_label = "Gravity" if is_gravity else "ÄKTA"
        elif "Knockout" in construct:
            base_label = "Knockout"
        else:
            base_label = "2Trig" if is_2trig else "1Trig"
        if is_solubility:
            label = f"{base_label} {time_label_from_test(test_name)}"
        else:
            label = base_label
    # Case 2: comparing 1 hour vs later timepoint for one construct
    elif len(protein_constructs) == 1 and has_1hr and has_later_timepoint:
        if "Knockout" in construct:
            trig_label = "Knockout"
        elif is_gravity:
            trig_label = "Gravity"
        else:
            trig_label = "2Trig" if is_2trig else "1Trig"
        label = f"{trig_label} {time_label_from_test(test_name)}"
    # Case 3: one replicate / one construct / one test
    elif len(protein_constructs) == 1 and len(tests) == 1:
        if "Knockout" in construct:
            base_label = "Knockout"
        elif is_gravity:
            base_label = "Gravity"
        else:
            base_label = "2Trig" if is_2trig else "1Trig"
        if is_solubility:
            label = f"{base_label} {time_label_from_test(test_name)}"
        else:
            label = base_label
    else:
        label = construct
    return f"{label} (n = {n})"

# This helper function creates the legend label for the inflection point line. 
def format_inflection_label(prefix, inflection_value):
    return f"{prefix}: pH = {inflection_value:.1f}"

# This helper function determines if the vWA of PA-UBA 1 Hour solubility data need to be labeld as behavior 1 (ideal) or behavior 2 (non-ideal).
def behavior_label_from_behavior_id(behavior_id):
    if behavior_id in ["W1", "P1"]:
        return "Behavior 1"
    elif behavior_id in ["W2", "P2"]:
        return "Behavior 2"
    return behavior_id

# This helper function assigns two unique colors to the data points and 4PL curve for the plots that show the different behaviors. It uses the base color assigned
# to behavior 1 as a reference for determining behavior 2's color.
def get_behavior_colors(base_color, behavior_ids):
    sorted_ids = sorted(behavior_ids)
    color_map = {}
    for i, behavior_id in enumerate(sorted_ids):
        factor = 1.0 - 0.35 * i
        if len(base_color) == 4:
            color = tuple(max(0, c * factor) for c in base_color[:3]) + (base_color[3],)
        else:
            color = tuple(max(0, c * factor) for c in base_color[:3])
        color_map[behavior_id] = color
    return color_map

# This helper function defines the y axis boundaries that the plot can show. We determined each boundary by looking at the plot and deciding how far zoomed in we needed to
# be in order to see the plot more clearly. This does visually exclude some points from the plot, but it doesn't affect how we calculate the 95% confidence interval or the
# 4PL curve (for plots that do have a 4PL curve).
def apply_custom_y_limits(protein_constructs, tests, plot_singly):
    key = (tuple(sorted(protein_constructs)), tuple(sorted(tests)))
    combined_plot_limits = {
        (tuple(sorted(("1TEL-PA-TNK1.UBA", "2Trig-1TEL-PA-TNK1.UBA"))), ("A400",)): {"ymax": 5.0},
        (tuple(sorted(("10xHis-1TEL-TV-vWA", "10xHis-1TEL-TV-vWA (Gravity)"))), ("A280_1hr",)): {"ymin": -1.0},
        (("1TEL-GG-TNK1.UBA",), ("A280_1hr", "A280_48-72hr")): {"ymin": -1.0},
        (("2Trig-1TEL-GG-TNK1.UBA",), ("A280_1hr", "A280_48-72hr")): {"ymin": -0.2},
        (("1TEL-PA-TNK1.UBA",), ("A280_1hr", "A280_48-72hr")): {"ymin": -0.4},
        (tuple(sorted(("2Trig-10xHis-1TEL-TV-vWA", "10xHis-1TEL-TV-vWA"))), ("A280_1hr",)): {"ymin": -0.8},
        (tuple(sorted(("2Trig-1TEL-GG-TNK1.UBA", "1TEL-GG-TNK1.UBA"))), ("A280_1hr",)): {"ymin": -0.4},
        (tuple(sorted(("2Trig-1TEL-PA-TNK1.UBA", "1TEL-PA-TNK1.UBA"))), ("A280_1hr",)): {"ymin": -0.2, "ymax": 1.7},
        (tuple(sorted(("2Trig-10xHis-1TEL-SR-TNK1.UBA", "10xHis-1TEL-SR-TNK1.UBA"))), ("A400",)): {"ymin": -0.2, "ymax": 1.8},}
    limits = None if plot_singly else combined_plot_limits.get(key)
    if limits is not None:
        plt.ylim(bottom=limits.get("ymin"), top=limits.get("ymax"))

def expand_tests_by_available_timepoints(proteins_and_tests):
    expanded = []
    for protein_construct, test in proteins_and_tests:
        if test != "A280_48-72hr":
            expanded.append((protein_construct, test))
            continue
        has_48hr = False
        has_72hr = False
        for bio_rep in protein_dict.get(protein_construct, {}):
            if len(protein_dict[protein_construct][bio_rep].get("A280_48hr", [])) > 0:
                has_48hr = True
            if len(protein_dict[protein_construct][bio_rep].get("A280_72hr", [])) > 0:
                has_72hr = True
        if has_48hr:
            expanded.append((protein_construct, "A280_48hr"))
        if has_72hr:
            expanded.append((protein_construct, "A280_72hr"))
    return expanded  

def statisticize(proteins_and_tests:list,plot_singly=False,plot=False,tailored=False, specific_date=None,specific_stock_id=None):
    data = {'protein construct':[],'date':[],'pH':[],'test':[],'data':[], 'is_starred':[], 'is_single_omitted':[], 'is_time_omitted':[], 'behavior_id':[]}
    proteins_and_tests = expand_tests_by_available_timepoints(proteins_and_tests)
    protein_constructs = list(dict.fromkeys(construct for construct, _ in proteins_and_tests))
    tests = list(dict.fromkeys(test for _, test in proteins_and_tests))
    is_alone_knockout_comparison = {
        "1TEL Alone", "2Trig-1TEL Alone", "1TEL Knockout"
    }.issubset(protein_constructs)
    is_time_comparison = (len(protein_constructs) == 1 and "A280_1hr" in tests and any(t in tests for t in ["A280_48hr", "A280_72hr"]))
    for protein_construct,test in proteins_and_tests:
        for bio_rep in protein_dict[protein_construct]:
            date = bio_rep.split("_")[0]
            if specific_date is not None and date != specific_date:
                continue
            stock_id = protein_dict[protein_construct][bio_rep].get('Stock ID', '').strip()
            if specific_stock_id is not None and stock_id != specific_stock_id:
                continue
            if is_time_comparison:
                has_1hr = len(protein_dict[protein_construct][bio_rep]["A280_1hr"]) != 0
                has_later_timepoint = (
                    len(protein_dict[protein_construct][bio_rep].get("A280_48hr", [])) != 0 or len(protein_dict[protein_construct][bio_rep].get("A280_72hr", [])) != 0)
                if not (has_1hr and has_later_timepoint):
                    continue
            if len(protein_dict[protein_construct][bio_rep][test])!=0:
                #if the biological replicate actually had this test done, aggregate its data into the dataframe replicate by replicate, conserving their order.
                for replicate in protein_dict[protein_construct][bio_rep][test]:
                    data['protein construct'].append(protein_construct)
                    data['date'].append(date)
                    pH = bio_rep.split("_")[-1]
                    data['pH'].append(pH)
                    data['test'].append(test)
                    behavior_id = protein_dict[protein_construct][bio_rep].get('Behavior ID', '').strip()
                    data['behavior_id'].append(behavior_id)
                    concentration = float(bio_rep.split("_")[1])
                    theoretical_max_concentration = concentration/5

                    # The "if isinstance" and it's associated "else" statement determines whether this data point needs to be omitted from the curve fitting. If it doesn't,
                    # then return False statements for certain labels. If it does need to be omitted, we return True statements that stick with the data point throughout
                    # the rest of the function. It also removes the symbol (*, ^, ~) used to signal that that data point needs to be omitted from curve fitting.
                    if isinstance(replicate, float):
                        numeric_replicate = replicate
                        is_starred = False
                        is_single_omitted = False
                        is_time_omitted = False
                    else:
                        replicate_str = replicate.strip()
                        is_starred = "*" in replicate_str
                        is_single_omitted = "^" in replicate_str
                        is_time_omitted = "~" in replicate_str
                        numeric_replicate = float(replicate_str.replace("*", "").replace("^", "").replace("~", ""))

                    # Normalized Precipitation Loss (solubility plots)
                    if test in ['A280_1hr', 'A280_48hr', 'A280_72hr']:
                        value = (theoretical_max_concentration - numeric_replicate) / theoretical_max_concentration
                    # Normalized Absorbance (turbidity plots)
                    else:
                        value = numeric_replicate * (5/concentration)
                    data['data'].append(value)
                    data['is_starred'].append(is_starred)
                    data['is_single_omitted'].append(is_single_omitted)
                    data['is_time_omitted'].append(is_time_omitted)

    #--------------------Now that the main filtering is done, it's time to make confidence intervals and fit curves for each set of construct:pH values---------------#

    constructs_by_pH_by_date = {}
    for replicate in range(len(data['protein construct'])):
        construct = data['protein construct'][replicate]
        test_name = data['test'][replicate]
        construct_key = (construct, test_name)
        pH = data['pH'][replicate]
        date = data['date'][replicate]
        measurement = data['data'][replicate]
        if construct_key not in constructs_by_pH_by_date:
            constructs_by_pH_by_date[construct_key] = {}
        is_starred = data['is_starred'][replicate]
        is_single_omitted = data['is_single_omitted'][replicate]
        is_time_omitted = data['is_time_omitted'][replicate]
        if pH not in constructs_by_pH_by_date[construct_key]:
            constructs_by_pH_by_date[construct_key][pH] = {"dates": {},"dates_starred": {},"dates_single_omitted": {},"dates_time_omitted": {},"dates_behavior_id": {}}
        if date not in constructs_by_pH_by_date[construct_key][pH]["dates"]:
            constructs_by_pH_by_date[construct_key][pH]["dates"][date] = []
            constructs_by_pH_by_date[construct_key][pH]["dates_starred"][date] = []
            constructs_by_pH_by_date[construct_key][pH]["dates_single_omitted"][date] = []
            constructs_by_pH_by_date[construct_key][pH]["dates_time_omitted"][date] = []
        constructs_by_pH_by_date[construct_key][pH]["dates"][date].append(measurement)
        constructs_by_pH_by_date[construct_key][pH]["dates_starred"][date].append(is_starred)
        constructs_by_pH_by_date[construct_key][pH]["dates_single_omitted"][date].append(is_single_omitted)
        constructs_by_pH_by_date[construct_key][pH]["dates_time_omitted"][date].append(is_time_omitted)
        constructs_by_pH_by_date[construct_key][pH]["dates_behavior_id"][date] = data['behavior_id'][replicate]

    biological_n_by_construct = {}

    for construct in constructs_by_pH_by_date:
        unique_dates = set()
        for pH in constructs_by_pH_by_date[construct]:
            for date in constructs_by_pH_by_date[construct][pH]["dates"]:
                unique_dates.add(date)
        biological_n_by_construct[construct] = len(unique_dates)

    biological_n_by_construct_and_behavior = {}
    for construct in constructs_by_pH_by_date:
        biological_n_by_construct_and_behavior[construct] = {}
        for pH in constructs_by_pH_by_date[construct]:
            for date in constructs_by_pH_by_date[construct][pH]["dates"]:
                behavior_id = constructs_by_pH_by_date[construct][pH]["dates_behavior_id"].get(date, "").strip()
                if behavior_id == "":
                    continue
                if behavior_id not in biological_n_by_construct_and_behavior[construct]:
                    biological_n_by_construct_and_behavior[construct][behavior_id] = set()
                biological_n_by_construct_and_behavior[construct][behavior_id].add(date)

    for construct in biological_n_by_construct_and_behavior:
        for behavior_id in biological_n_by_construct_and_behavior[construct]:
            biological_n_by_construct_and_behavior[construct][behavior_id] = len(biological_n_by_construct_and_behavior[construct][behavior_id])

    # Compute pointwise 95% confidence intervals from the biological replicate dates.
    for construct in constructs_by_pH_by_date:
        for pH in constructs_by_pH_by_date[construct]:
            pH_dates = constructs_by_pH_by_date[construct][pH]["dates"]
            pH_starred = constructs_by_pH_by_date[construct][pH]["dates_starred"]
            pH_single_omitted = constructs_by_pH_by_date[construct][pH]["dates_single_omitted"]
            pH_time_omitted = constructs_by_pH_by_date[construct][pH]["dates_time_omitted"]

            # Each date contributes one observation, regardless of how many drops were measured.
            date_means_for_display = []
            date_means_for_fit = []

            # This for loop goes through each pH for one construct and "attaches" the labels we use to determine if a measurement should be omitted from curve fitting
            # and if it should have an asterisk next to it in the plot signifying that it was omitted from curve fitting.
            for date_key in pH_dates:
                drops = pH_dates[date_key]
                star_flags = pH_starred[date_key]
                single_omit_flags = pH_single_omitted[date_key]
                time_omit_flags = pH_time_omitted[date_key]
                paired = list(zip(drops, star_flags, single_omit_flags, time_omit_flags))
                fit_paired = paired

                # d = drop, s = star flag (omitted from combined plots), so = just omit from single plot, to = just omit if doing a 1 Hour vs 48 Hour plot.
                # Remove starred points (if tailored) (csv symbol: *)
                if tailored:
                    fit_paired = [item for item in fit_paired if not item[1]]
                # Remove single-plot omissions (csv symbol: ^)
                if plot_singly:
                    fit_paired = [item for item in fit_paired if not item[2]]
                # Remove 1 Hour vs 48 Hour omissions (csv symbol: ~)
                if is_time_comparison:
                    fit_paired = [item for item in fit_paired if not item[3]]
                display_drops = [item[0] for item in paired]
                fit_drops = [item[0] for item in fit_paired]
                if len(display_drops) > 0:
                    date_means_for_display.append(np.mean(display_drops))
                if len(fit_drops) > 0:
                    date_means_for_fit.append(np.mean(fit_drops))

            # Create variables that will be displayed on the plot and used for calculating the 95% confidence interval.
            if date_means_for_display:
                mean_display = np.mean(date_means_for_display)
                constructs_by_pH_by_date[construct][pH]["mean_display"] = mean_display
                if len(date_means_for_display) > 1:
                    sd_display = np.std(date_means_for_display, ddof=1)
                    se_display = sd_display / np.sqrt(len(date_means_for_display))
                    t_display = t.ppf(0.975, df=len(date_means_for_display) - 1)
                    error_display = se_display * t_display
                    constructs_by_pH_by_date[construct][pH]["CI_lower_display"] = mean_display - error_display
                    constructs_by_pH_by_date[construct][pH]["CI_upper_display"] = mean_display + error_display

            # Calculate the mean of the data we use to do the curve fitting.
            if date_means_for_fit:
                constructs_by_pH_by_date[construct][pH]["mean"] = np.mean(date_means_for_fit)
    
    # This function returns the 4PL equation that we use to plot the 4PL curve onto the plot. We get the parameters from the "curve_fit" command below.
    def pH_to_absorbance_model_4pl(pH,upper_asymptote,Hill_slope,inflection_point,lower_asymptote):
        return lower_asymptote + (upper_asymptote - lower_asymptote) / (1 + (pH / inflection_point)**Hill_slope)

    # This helper function categorizes datasets in combined plots. Since all combined plots are comparing one group to another (1Trig vs 2Trig, 1 Hour vs 48 Hour, Gravity vs ÄKTA).
    # Group 1 will be placed higher than group 2 in the legend.
    def classify_legend_group(construct, test_name, protein_constructs, tests):
        is_gravity_comparison = (len(protein_constructs) == 2 and len(tests) == 1 and any("Gravity" in c for c in protein_constructs))
        is_time_comparison = (len(protein_constructs) == 1 and "A280_1hr" in tests and any(t in tests for t in ["A280_48hr", "A280_72hr"]))
        if is_gravity_comparison:
            if "Gravity" in construct:
                return "group2"
            return "group1"
        if is_time_comparison:
            if test_name in ["A280_48hr", "A280_72hr"]:
                return "group2"
            return "group1"
        if "2Trig" in construct:
            return "group2"
        return "group1"
    
    if plot==True:
        plt.figure(figsize=(8.75,6))
        handles_and_labels = {}
        inflection_handles_and_labels = {}
        legend_groups = {"group1_main": [],"group1_inflection": [],"group1_ci": [],"group2_main": [],"group2_inflection": [],"group2_ci": [],"other": []}
        combined_any_omitted_points = False

        def change_colors():
            keys_1trig = []
            keys_2trig = []
            keys_knockout = []
            keys_48_hour = []
            keys_2Trig_48_hour = []
            keys_Gravity = []
            keys_Gravity_48_hour = []
            for construct_key in constructs_by_pH_by_date:
                construct, test_name = construct_key
                if "Knockout" in construct:
                    keys_knockout.append(construct_key)
                elif "Gravity" in construct and test_name in ["A280_48hr", "A280_72hr"]:
                    keys_Gravity_48_hour.append(construct_key)
                elif "Gravity" in construct:
                    keys_Gravity.append(construct_key)
                elif test_name in ["A280_48hr", "A280_72hr"] and not '2Trig' in construct:
                    keys_48_hour.append(construct_key)
                elif test_name in ["A280_48hr", "A280_72hr"] and "2Trig" in construct:
                    keys_2Trig_48_hour.append(construct_key)
                elif "2Trig" in construct:
                    keys_2trig.append(construct_key)
                else:
                    keys_1trig.append(construct_key)

            def shade_list(cmap, n, lo=0.25, hi=0.9):
                if n<= 1:
                    return [cmap((lo + hi) / 2)]
                return [cmap(x) for x in np.linspace(lo, hi, n)]

            # Lists of color shades for the 95% confidence intervals.
            blue_shades = shade_list(plt.cm.Blues, len(keys_1trig))
            orange_shades = shade_list(plt.cm.Oranges, len(keys_2trig))
            green_shades = shade_list(plt.cm.Greens, len(keys_knockout))
            crimson_shades = shade_list(plt.cm.PuRd, len(keys_48_hour))
            purple_shades = shade_list(plt.cm.Purples, len(keys_2Trig_48_hour))
            lime_shades = shade_list(plt.cm.YlGn, len(keys_Gravity))
            brown_shades = shade_list(plt.cm.copper, len(keys_Gravity_48_hour))

            # Determines colors for the scatter plots.
            color_map = {}
            for construct, color in zip(sorted(keys_1trig), blue_shades):
                color_map[construct] = color
            for construct, color in zip(sorted(keys_2trig), orange_shades):
                color_map[construct] = color
            for construct, color in zip(sorted(keys_knockout), green_shades):
                color_map[construct] = color
            for construct, color in zip(sorted(keys_48_hour), crimson_shades):
                color_map[construct] = color
            for construct, color in zip(sorted(keys_2Trig_48_hour), purple_shades):
                color_map[construct] = color
            for construct, color in zip(sorted(keys_Gravity), lime_shades):
                color_map[construct] = color
            for construct, color in zip(sorted(keys_Gravity_48_hour), brown_shades):
                color_map[construct] = color
            return color_map
        color_map = change_colors()

        inflection_points_1trig_1hr = []
        inflection_points_1trig_48hr = []
        inflection_points_2trig_1hr = []
        inflection_points_2trig_48hr = []
        inflection_points_knockout_1hr = []
        inflection_points_knockout_48hr = []
        inflection_points_gravity_1hr = []
        inflection_points_gravity_48hr = []

        for construct_key in constructs_by_pH_by_date:
            any_starred_points = False
            construct, test_name = construct_key
            is_special_behavior_plot = (
                test_name == "A280_1hr"
                and len(protein_constructs) == 2
                and len(tests) == 1
                and construct in {"10xHis-1TEL-TV-vWA", "1TEL-PA-TNK1.UBA"}
            )
            pH_values_fit = []
            means_fit = []
            pH_values_display = []
            ci_lower_bounds = []
            ci_upper_bounds = []
            construct_pHs =[]
            construct_data = []
            construct_starred = []
            construct_single_omitted = []
            construct_behavior_ids = []
            construct_time_omitted = []
            special_curve_handles = {}
            special_scatter_handles = {}
            # Sort pHs numerically.
            sorted_pHs = sorted(constructs_by_pH_by_date[construct_key].keys(), key=float)
            behavior_groups = {}
            for pH in sorted_pHs:
                pH_entry = constructs_by_pH_by_date[construct_key][pH]
                if is_special_behavior_plot:
                    for date in pH_entry["dates"]:
                        behavior_id = pH_entry["dates_behavior_id"].get(date, "").strip()
                        if behavior_id == "":
                            continue
                        if behavior_id not in behavior_groups:
                            behavior_groups[behavior_id] = {
                                'pH_values': [],
                                'means': [],
                                'ci_pH_values': [],
                                'ci_lower': [],
                                'ci_upper': [],
                            }
                # Always keep raw data points for plotting.
                for date, measurements in pH_entry["dates"].items():
                    star_flags = pH_entry["dates_starred"][date]
                    single_omit_flags = pH_entry["dates_single_omitted"][date]
                    time_omit_flags = pH_entry["dates_time_omitted"][date]
                    behavior_id = pH_entry["dates_behavior_id"].get(date, "").strip()
                    for measurement, is_starred, is_single_omitted, is_time_omitted in zip(measurements, star_flags, single_omit_flags, time_omit_flags):
                        construct_pHs.append(float(pH))
                        construct_data.append(measurement)
                        construct_starred.append(is_starred)
                        construct_single_omitted.append(is_single_omitted)
                        construct_behavior_ids.append(behavior_id)
                        construct_time_omitted.append(is_time_omitted)
                # If we weren't able to compute statistics for data at a certain pH then we don't plot the confidence interval bands or any of the scatter points.
                if "mean" not in pH_entry and "mean_display" not in pH_entry:
                    continue
                # The following if and else statements compute statistics for datasets with unique behaviors and add those data to empty lists. We also add any other data
                # to other lists.
                if is_special_behavior_plot:
                    subgroup_measurements_by_date = {}
                    for date, measurements in pH_entry["dates"].items():
                        behavior_id = pH_entry["dates_behavior_id"].get(date, "").strip()
                        if behavior_id == "":
                            continue
                        if behavior_id not in subgroup_measurements_by_date:
                            subgroup_measurements_by_date[behavior_id] = {}
                        star_flags = pH_entry["dates_starred"][date]
                        single_omit_flags = pH_entry["dates_single_omitted"][date]
                        time_omit_flags = pH_entry["dates_time_omitted"][date]
                        paired = list(zip(measurements, star_flags, single_omit_flags, time_omit_flags))
                        if tailored:
                            paired = [item for item in paired if not item[1]]
                        if plot_singly:
                            paired = [item for item in paired if not item[2]]
                        if is_time_comparison:
                            paired = [item for item in paired if not item[3]]
                        if paired:
                            subgroup_measurements_by_date[behavior_id][date] = np.mean(
                                [item[0] for item in paired])
                    for behavior_id, date_means in subgroup_measurements_by_date.items():
                        if not date_means:
                            continue
                        vals = np.array(list(date_means.values()), dtype=float)
                        mean_val = np.mean(vals)
                        behavior_groups[behavior_id]['pH_values'].append(float(pH))
                        behavior_groups[behavior_id]['means'].append(mean_val)
                        if len(vals) > 1:
                            sd_val = np.std(vals, ddof=1)
                            se_val = sd_val / np.sqrt(len(vals))
                            t_val = t.ppf(0.975, df=len(vals) - 1)
                            err = se_val * t_val
                            behavior_groups[behavior_id]['ci_pH_values'].append(float(pH))
                            behavior_groups[behavior_id]['ci_lower'].append(mean_val - err)
                            behavior_groups[behavior_id]['ci_upper'].append(mean_val + err)
                else:
                    # Add CI-band values whenever display stats exist
                    if "CI_lower_display" in pH_entry:
                        pH_values_display.append(float(pH))
                        ci_lower_bounds.append(pH_entry["CI_lower_display"])
                        ci_upper_bounds.append(pH_entry["CI_upper_display"])

                    # Add fit values only when fit stats exist
                    if "mean" in pH_entry:
                        pH_values_fit.append(float(pH))
                        means_fit.append(pH_entry["mean"])

        #--------------------The rest of the code at this point until "if __name__ == '__main__':" deals with curve fitting and plotting---------------------#

            if is_special_behavior_plot:
                behavior_ids_present = sorted([bid for bid in behavior_groups.keys() if bid != ""])
                behavior_colors = get_behavior_colors(color_map[construct_key], behavior_ids_present)

            try:
                if is_special_behavior_plot:
                    for behavior_id, group_data in behavior_groups.items():
                        if len(group_data['pH_values']) < 3:
                            continue
                        group_pH = group_data['pH_values']
                        group_means = group_data['means']
                        group_ci_pH = group_data['ci_pH_values']
                        group_ci_lower = group_data['ci_lower']
                        group_ci_upper = group_data['ci_upper']
                        pH_linspace = np.linspace(min(group_pH), max(group_pH), 400)
                        initial_guesses = [max(group_means), -4, 5.5, min(group_means)]
                        # This is where we fit a curve to the data using scipy's curve_fit function. Scipy uses non-linear least squares to fit a function, f, to data. 
                        best_fit_parameters, _ = curve_fit(pH_to_absorbance_model_4pl,group_pH,group_means,p0=initial_guesses)
                        inflection_value = round(best_fit_parameters[2],1)
                        if inflection_value <= 4.0:
                            continue
                        curve_handle_list = plt.plot(pH_linspace,pH_to_absorbance_model_4pl(pH_linspace, *best_fit_parameters),color=behavior_colors[behavior_id],linestyle="-")
                        line_handle = curve_handle_list[0]
                        if group_ci_pH:
                            ci_handle = plt.fill_between(
                                group_ci_pH, group_ci_lower, group_ci_upper,
                                color=behavior_colors[behavior_id], alpha=0.12)
                            handles_and_labels[ci_handle] = f"{behavior_label_from_behavior_id(behavior_id)} 95% CI"
                            legend_groups["group1_ci"].append(ci_handle)
                        inflection_handle = plt.axvline(inflection_value,color=behavior_colors[behavior_id],linestyle=":")
                        inflection_handles_and_labels[inflection_handle] = (
                            f"{behavior_label_from_behavior_id(behavior_id)} Inflection: "
                            f"{format_inflection_label('', inflection_value).replace(': ', '').strip()}")
                        legend_groups["group1_inflection"].append(inflection_handle)
                        special_curve_handles[behavior_id] = line_handle
                else:
                    pH_linspace = np.linspace(min(pH_values_fit), max(pH_values_fit), 400)
                    initial_guesses = [max(means_fit), -4, 5.5, min(means_fit)]
                    best_fit_parameters, _ = curve_fit(pH_to_absorbance_model_4pl, pH_values_fit, means_fit, p0=initial_guesses)
                    inflection_value = round(best_fit_parameters[2],1)
                    if inflection_value <= 4.0:
                        raise ValueError(f"Calculated inflection point ({inflection_value:.1f}) is <= 4.0")
                    if "Knockout" in construct:
                        if test_name in ["A280_48hr", "A280_72hr"]:
                            inflection_points_knockout_48hr.append(inflection_value)
                        else:
                            inflection_points_knockout_1hr.append(inflection_value)
                    elif "Gravity" in construct:
                        if test_name in ["A280_48hr", "A280_72hr"]:
                            inflection_points_gravity_48hr.append(inflection_value)
                        else:
                            inflection_points_gravity_1hr.append(inflection_value)
                    elif "2Trig" in construct:
                        if test_name in ["A280_48hr", "A280_72hr"]:
                            inflection_points_2trig_48hr.append(inflection_value)
                        else:
                            inflection_points_2trig_1hr.append(inflection_value)
                    else:
                        if test_name in ["A280_48hr", "A280_72hr"]:
                            inflection_points_1trig_48hr.append(inflection_value)
                        else:
                            inflection_points_1trig_1hr.append(inflection_value)
                    curve_handle_list = plt.plot(pH_linspace,pH_to_absorbance_model_4pl(pH_linspace, *best_fit_parameters),color=color_map[construct_key])
                    line_handle = curve_handle_list[0]
                    if pH_values_display:
                        if "Gravity" in construct and test_name in ["A280_48hr", "A280_72hr"]:
                            ci_handle = plt.fill_between(pH_values_display, ci_lower_bounds, ci_upper_bounds, color='brown', alpha=0.20)
                            ci_label = f"Gravity {time_label_from_test(test_name)} 95% CI"
                        elif "Gravity" in construct and test_name == 'A280_1hr':
                            ci_handle = plt.fill_between(pH_values_display, ci_lower_bounds, ci_upper_bounds, color='lime', alpha=0.20)
                            ci_label = "Gravity 1 Hour 95% CI"
                        elif "Gravity" in construct:
                            ci_handle = plt.fill_between(pH_values_display, ci_lower_bounds, ci_upper_bounds, color='lime', alpha=0.20)
                            ci_label = "Gravity 95% CI"
                        elif "Knockout" in construct:
                            ci_handle = plt.fill_between(pH_values_display, ci_lower_bounds, ci_upper_bounds, color=color_map[construct_key], alpha=0.20)
                            time_label = time_label_from_test(test_name)
                            ci_label = f"Knockout {time_label + ' ' if time_label else ''}95% CI"
                        elif test_name in ["A280_48hr", "A280_72hr"] and "2Trig" not in construct:
                            ci_handle = plt.fill_between(pH_values_display, ci_lower_bounds, ci_upper_bounds, color='pink', alpha=0.20)
                            ci_label = f"1Trig {time_label_from_test(test_name)} 95% CI"
                        elif test_name in ["A280_48hr", "A280_72hr"] and "2Trig" in construct:
                            ci_handle = plt.fill_between(pH_values_display, ci_lower_bounds, ci_upper_bounds, color='mediumpurple', alpha=0.20)
                            ci_label = f"2Trig {time_label_from_test(test_name)} 95% CI"
                        elif "2Trig" in construct and test_name == 'A280_1hr':
                            ci_handle = plt.fill_between(pH_values_display, ci_lower_bounds, ci_upper_bounds, color='tomato', alpha=0.20)
                            ci_label = "2Trig 1 Hour 95% CI"
                        elif "2Trig" in construct:
                            ci_handle = plt.fill_between(pH_values_display, ci_lower_bounds, ci_upper_bounds, color='tomato', alpha=0.20)
                            ci_label = "2Trig 95% CI"
                        elif test_name == 'A280_1hr' and "2Trig" not in construct:
                            ci_handle = plt.fill_between(pH_values_display, ci_lower_bounds, ci_upper_bounds, color='cyan', alpha=0.20)
                            ci_label = "1Trig 1 Hour 95% CI"
                        else:
                            ci_handle = plt.fill_between(pH_values_display, ci_lower_bounds, ci_upper_bounds, color='cyan', alpha=0.20)
                            ci_label = "1Trig 95% CI"
                        handles_and_labels[ci_handle] = ci_label
                        group_prefix = classify_legend_group(construct, test_name, protein_constructs, tests)
                        legend_groups[f"{group_prefix}_ci"].append(ci_handle)

                # plot raw scatter points
                jitter_strength = 0.08
                jittered_pHs = np.array(construct_pHs) + np.random.uniform(-jitter_strength, jitter_strength,len(construct_pHs))
                if is_special_behavior_plot:
                    special_scatter_handles = {}
                    behavior_ids_present = sorted(set([bid for bid in construct_behavior_ids if bid != ""]))
                    behavior_colors = get_behavior_colors(color_map[construct_key], behavior_ids_present)

                    available_markers = ["^", "s", "o", "D", "v", "P", "X"]
                    marker_map = {bid: available_markers[i % len(available_markers)] for i, bid in enumerate(behavior_ids_present)}

                    for behavior_id, marker_style in marker_map.items():
                        xs = [x for x, gid in zip(jittered_pHs, construct_behavior_ids) if gid == behavior_id]
                        ys = [y for y, gid in zip(construct_data, construct_behavior_ids) if gid == behavior_id]
                        if len(xs) > 0:
                            scatter_handle = plt.scatter(xs, ys, marker=marker_style, color=behavior_colors[behavior_id], alpha=0.2)
                            special_scatter_handles[behavior_id] = scatter_handle

                    for behavior_id in special_curve_handles:
                        line_handle = special_curve_handles[behavior_id]
                        scatter_handle = special_scatter_handles.get(behavior_id)
                        behavior_label = behavior_label_from_behavior_id(behavior_id)
                        subgroup_n = biological_n_by_construct_and_behavior.get(construct_key, {}).get(behavior_id, 0)
                        base_label = format_legend_label(construct, test_name, protein_constructs, tests, subgroup_n)

                        if scatter_handle is not None:
                            pair_handle = (line_handle, scatter_handle)
                            handles_and_labels[pair_handle] = f"{base_label}, {behavior_label}"
                            legend_groups["group1_main"].append(pair_handle)
                        else:
                            handles_and_labels[line_handle] = f"{base_label}, {behavior_label}"
                            legend_groups["group1_main"].append(line_handle)
                elif "Knockout" in construct:
                    scatter_handle = plt.scatter(jittered_pHs, construct_data, marker="s", color=color_map[construct_key], alpha=0.3)
                elif "2Trig" in construct and test_name == 'A280_48-72hr':
                    scatter_handle = plt.scatter(jittered_pHs, construct_data, marker="^",color=color_map[construct_key], alpha=0.3)
                elif "2Trig" in construct:
                    scatter_handle = plt.scatter(jittered_pHs, construct_data, marker="^",color=color_map[construct_key], alpha=0.2)
                else:
                    scatter_handle = plt.scatter(jittered_pHs, construct_data, color=color_map[construct_key],alpha=0.2)
                # Plot asterisks for omitted (starred) data
                omitted_flags = [
                    s or (plot_singly and so) or (is_time_comparison and to)
                    for s, so, to in zip(construct_starred, construct_single_omitted, construct_time_omitted)]
                omitted_jittered_pHs = [x for x, omit in zip(jittered_pHs, omitted_flags) if omit]
                omitted_data = [y for y, omit in zip(construct_data, omitted_flags) if omit]
                if len(omitted_data) > 0:
                    for x, y in zip(omitted_jittered_pHs, omitted_data):
                        plt.annotate("*",(x, y),xytext=(0.005, 0.00005),textcoords="offset points",ha="left",va="bottom",color="black",fontsize=11,fontweight="light",zorder=6)
                    any_starred_points = True
                    combined_any_omitted_points = True
                if not is_special_behavior_plot:
                    legend_label = format_legend_label(construct,test_name,protein_constructs,tests,biological_n_by_construct[construct_key])
                    pair_handle = (line_handle, scatter_handle)
                    handles_and_labels[pair_handle] = legend_label
                    group_prefix = classify_legend_group(construct, test_name, protein_constructs, tests)
                    legend_groups[f"{group_prefix}_main"].append(pair_handle)

                if plot_singly:
                    if "Gravity" in construct and test_name in ["A280_48hr", "A280_72hr"]:
                        inflection_color = 'brown'
                    elif "Gravity" in construct:
                        inflection_color = "green"
                    elif "Knockout" in construct:
                        inflection_color = "forestgreen"
                    elif test_name in ["A280_48hr", "A280_72hr"] and "2Trig" in construct:
                        inflection_color = "indigo"
                    elif test_name in ["A280_48hr", "A280_72hr"] and "2Trig" not in construct:
                        inflection_color = "purple"
                    elif "2Trig" in construct:
                        inflection_color = "red"
                    else:
                        inflection_color = "blue"
                    inflection_handle = plt.axvline(inflection_value,color=inflection_color,linestyle="--")
                    handles_and_labels[inflection_handle] = format_inflection_label("Inflection Point",inflection_value)

            except Exception as e:
                print(f'The curve fit for {construct} failed due to: {e}. Plotting points instead of fitted curve...')
                construct_pHs = []
                construct_data = []
                construct_starred = []
                construct_single_omitted = []
                construct_time_omitted = []
                sorted_pHs = sorted(constructs_by_pH_by_date[construct_key].keys(), key=float)
                for pH in sorted_pHs:
                    pH_entry = constructs_by_pH_by_date[construct_key][pH]
                    for date, measurements in pH_entry["dates"].items():
                        star_flags = pH_entry["dates_starred"][date]
                        single_omit_flags = pH_entry["dates_single_omitted"][date]
                        time_omit_flags = pH_entry["dates_time_omitted"][date]
                        for measurement, is_starred, is_single_omitted, is_time_omitted in zip(measurements, star_flags, single_omit_flags,time_omit_flags):
                            construct_pHs.append(float(pH))
                            construct_data.append(measurement)
                            construct_starred.append(is_starred)
                            construct_single_omitted.append(is_single_omitted)
                            construct_time_omitted.append(is_time_omitted)
                ci_handle = None
                if pH_values_display:
                    if "Gravity" in construct:
                        ci_handle = plt.fill_between(pH_values_display, ci_lower_bounds, ci_upper_bounds, color='lime', alpha=0.20)
                    elif "Knockout" in construct:
                        ci_handle = plt.fill_between(pH_values_display, ci_lower_bounds, ci_upper_bounds, color=color_map[construct_key], alpha=0.20)
                    elif test_name in ["A280_48hr", "A280_72hr"] and "2Trig" not in construct:
                        ci_handle = plt.fill_between(pH_values_display, ci_lower_bounds, ci_upper_bounds, color='pink', alpha=0.20)
                    elif test_name in ["A280_48hr", "A280_72hr"] and "2Trig" in construct:
                        ci_handle = plt.fill_between(pH_values_display, ci_lower_bounds, ci_upper_bounds, color='mediumpurple', alpha=0.20)
                    elif "2Trig" in construct:
                        ci_handle = plt.fill_between(pH_values_display, ci_lower_bounds, ci_upper_bounds, color='tomato', alpha=0.20)
                    else:
                        ci_handle = plt.fill_between(pH_values_display, ci_lower_bounds, ci_upper_bounds, color='cyan', alpha=0.20)
                jitter_strength = 0.08
                jittered_pHs = np.array(construct_pHs) + np.random.uniform(-jitter_strength, jitter_strength, len(construct_pHs))
                if "Knockout" in construct:
                    scatter_handle = plt.scatter(jittered_pHs, construct_data, marker="s", color=color_map[construct_key], alpha=0.3)
                elif "2Trig" in construct and test_name == 'A280_48-72hr':
                    scatter_handle = plt.scatter(jittered_pHs,construct_data,marker="^", color=color_map[construct_key], alpha=0.3)
                elif "2Trig" in construct:
                    scatter_handle = plt.scatter(jittered_pHs,construct_data,marker="^", color=color_map[construct_key], alpha=0.2)
                else:
                    scatter_handle = plt.scatter(jittered_pHs,construct_data, color=color_map[construct_key], alpha=0.2)
                omitted_flags = [
                    s or (plot_singly and so) or (is_time_comparison and to)
                    for s, so, to in zip(construct_starred, construct_single_omitted, construct_time_omitted)]
                omitted_jittered_pHs = [x for x, omit in zip(jittered_pHs, omitted_flags) if omit]
                omitted_data = [y for y, omit in zip(construct_data, omitted_flags) if omit]
                if len(omitted_data) > 0:
                    for x, y in zip(omitted_jittered_pHs, omitted_data):
                        plt.annotate("*",(x, y),xytext=(0.005, 0.00005),textcoords="offset points",ha="left",va="bottom",color="black",fontsize=11,fontweight="light",zorder=6)
                    any_starred_points = True
                    combined_any_omitted_points = True
                legend_label = format_legend_label(construct,test_name,protein_constructs,tests,biological_n_by_construct[construct_key])
                handles_and_labels[scatter_handle] = legend_label
                if ci_handle is not None:
                    handles_and_labels[ci_handle] = f"{legend_label.split(' (n')[0]} 95% CI"
                group_prefix = classify_legend_group(construct, test_name, protein_constructs, tests)
                if ci_handle is not None:
                    legend_groups[f"{group_prefix}_ci"].append(ci_handle)
                legend_groups[f"{group_prefix}_main"].append(scatter_handle)
            if plot_singly:
                if test_name == "A400":
                    y_axis = "Normalized Absorbance Units"
                    assay = "Turbidity"
                else:
                    y_axis = "Normalized Sedimentation"
                    assay = "Sedimentation"
                replicate_dates = sorted({d
                    for pH in constructs_by_pH_by_date[construct_key]
                    for d in constructs_by_pH_by_date[construct_key][pH]["dates"]})
                date_text = replicate_dates[0] if len(replicate_dates) == 1 else ", ".join(replicate_dates)
                plt.title(f"{format_chart_title(construct,len(protein_constructs))}, {assay}\nReplicate Date: {date_text}",fontsize=14)
                plt.ylabel(y_axis)
                plt.xlabel("pH")
                if any_starred_points:
                    asterisk_handle = plt.scatter([], [], marker='$*$', color='black', s=40)
                    handles_and_labels[asterisk_handle] = "Omitted from curve fitting"
                ordered_handles = []
                # 1) main data entry first (line+dot pair, or scatter-only fallback)
                for h in handles_and_labels:
                    label = handles_and_labels[h]
                    if "(n =" in label:
                        ordered_handles.append(h)
                # 2) inflection line(s) next
                for h in handles_and_labels:
                    label = handles_and_labels[h]
                    if "Inflection" in label:
                        ordered_handles.append(h)
                # 3) CI band next
                for h in handles_and_labels:
                    label = handles_and_labels[h]
                    if "95% CI" in label:
                        ordered_handles.append(h)
                # 4) omitted-points key last
                for h in handles_and_labels:
                    label = handles_and_labels[h]
                    if "Omitted from curve fitting" in label:
                        ordered_handles.append(h)
                # remove duplicates while preserving order
                seen = set()
                handles = []
                for h in ordered_handles:
                    if h not in seen:
                        handles.append(h)
                        seen.add(h)
                labels = [handles_and_labels[h] for h in handles]
                apply_custom_y_limits(protein_constructs, tests, plot_singly=True)
                ax = plt.gca()
                ax.yaxis.set_major_formatter(mticker.FormatStrFormatter('%.2f'))
                plt.legend(handles,labels,handler_map={tuple: HandlerTuple(ndivide=2, pad=1)},loc="upper right")
                safe_date_text = date_text.replace("/", "-").replace(", ", "_").replace(" ", "_")
                plt.savefig(f"{BASE_DIR}/figures_with_CIs/{construct.replace('/', '_')}_{test_name}_{safe_date_text}_single.png",dpi=300,bbox_inches="tight")
                plt.show()
                plt.close()
                handles_and_labels = {}

        if not plot_singly:
            is_gravity_comparison = (len(protein_constructs) == 2 and len(tests) == 1 and any("Gravity" in construct for construct in protein_constructs))
            is_time_comparison = (len(protein_constructs) == 1 and "A280_1hr" in tests and any(t in tests for t in ["A280_48hr", "A280_72hr"]))
            if is_time_comparison:
                if "Knockout" in protein_constructs[0]:
                    inflection_groups = [
                        (inflection_points_knockout_1hr, "forestgreen", "1 Hour Inflection"),
                        (inflection_points_knockout_48hr, "darkgreen", f"{time_label_from_test(test_name)} Inflection"),
                    ]
                    for inflection_points, color, label in inflection_groups:
                        if inflection_points:
                            mean_inflection = np.mean(inflection_points)
                            inflection_handle = plt.axvline(
                                mean_inflection, color=color, linestyle="--")
                            inflection_handles_and_labels[inflection_handle] = format_inflection_label(
                                label, mean_inflection)
                            legend_groups["group1_inflection" if "1 Hour" in label else "group2_inflection"].append(inflection_handle)
                elif "2Trig" in protein_constructs[0]:
                    if len(inflection_points_2trig_1hr) >= 1:
                        mean_1hr = np.mean(inflection_points_2trig_1hr)
                        inflection_1hr = plt.axvline(mean_1hr, color="darkorange", linestyle="--")
                        inflection_handles_and_labels[inflection_1hr] = format_inflection_label("1 Hour Inflection",mean_1hr)
                        legend_groups["group1_inflection"].append(inflection_1hr)

                    if len(inflection_points_2trig_48hr) >= 1:
                        mean_48hr = np.mean(inflection_points_2trig_48hr)
                        inflection_48hr = plt.axvline(mean_48hr, color="indigo", linestyle="--")
                        inflection_handles_and_labels[inflection_48hr] = format_inflection_label(f"{time_label_from_test(test_name)} Inflection",mean_48hr)
                        legend_groups["group2_inflection"].append(inflection_48hr)

                elif "Gravity" in protein_constructs[0]:
                    if len(inflection_points_gravity_1hr) >= 1:
                        mean_1hr = np.mean(inflection_points_gravity_1hr)
                        inflection_1hr = plt.axvline(mean_1hr, color="green", linestyle="--")
                        inflection_handles_and_labels[inflection_1hr] = format_inflection_label("1 Hour Inflection", mean_1hr)
                        legend_groups["group1_inflection"].append(inflection_1hr)

                    if len(inflection_points_gravity_48hr) >= 1:
                        mean_48hr = np.mean(inflection_points_gravity_48hr)
                        inflection_48hr = plt.axvline(mean_48hr, color="brown", linestyle="--")
                        inflection_handles_and_labels[inflection_48hr] = format_inflection_label(f"{time_label_from_test(test_name)} Inflection", mean_48hr)
                        legend_groups["group2_inflection"].append(inflection_48hr)   
                
                else:
                    if len(inflection_points_1trig_1hr) >= 1:
                        mean_1hr = np.mean(inflection_points_1trig_1hr)
                        inflection_1hr = plt.axvline(mean_1hr, color="blue", linestyle="--")
                        inflection_handles_and_labels[inflection_1hr] = format_inflection_label("1 Hour Inflection",mean_1hr)
                        legend_groups["group1_inflection"].append(inflection_1hr)

                    if len(inflection_points_1trig_48hr) >= 1:
                        mean_48hr = np.mean(inflection_points_1trig_48hr)
                        inflection_48hr = plt.axvline(mean_48hr, color="purple", linestyle="--")
                        inflection_handles_and_labels[inflection_48hr] = format_inflection_label(f"{time_label_from_test(test_name)} Inflection",mean_48hr)
                        legend_groups["group2_inflection"].append(inflection_48hr)
            elif is_gravity_comparison:
                if "A400" in tests:
                    nongravity_points = inflection_points_1trig_1hr
                    nongravity_color = "blue"
                    nongravity_label = "ÄKTA Inflection"

                    gravity_points = inflection_points_gravity_1hr
                    gravity_color = "green"
                    gravity_label = "Gravity Inflection"

                elif any(t in tests for t in ["A280_48hr", "A280_72hr"]):
                    nongravity_points = inflection_points_1trig_48hr
                    nongravity_color = "purple"
                    nongravity_label = f"ÄKTA {time_label_from_test(test_name)} Inflection"

                    gravity_points = inflection_points_gravity_48hr
                    gravity_color = "brown"
                    gravity_label = f"Gravity {time_label_from_test(test_name)} Inflection"

                elif "A280_1hr" in tests:
                    nongravity_points = inflection_points_1trig_1hr
                    nongravity_color = "blue"
                    nongravity_label = "ÄKTA 1 Hour Inflection"

                    gravity_points = inflection_points_gravity_1hr
                    gravity_color = "green"
                    gravity_label = "Gravity 1 Hour Inflection"

                if len(nongravity_points) >= 1:
                    mean_ng = np.mean(nongravity_points)
                    inflection_ng = plt.axvline(mean_ng, color=nongravity_color, linestyle="--")
                    inflection_handles_and_labels[inflection_ng] = format_inflection_label(nongravity_label, mean_ng)
                    legend_groups["group1_inflection"].append(inflection_ng)

                if len(gravity_points) >= 1:
                    mean_g = np.mean(gravity_points)
                    inflection_g = plt.axvline(mean_g, color=gravity_color, linestyle="--")
                    inflection_handles_and_labels[inflection_g] = format_inflection_label(gravity_label, mean_g)
                    legend_groups["group2_inflection"].append(inflection_g)
            else:
                # standard 1Trig vs 2Trig comparison
                if "A280_48-72hr" in tests:
                    inflection_points_1trig = inflection_points_1trig_48hr
                    inflection_points_2trig = inflection_points_2trig_48hr
                    color_1trig = "purple"
                    color_2trig = "indigo"
                else:
                    inflection_points_1trig = inflection_points_1trig_1hr
                    inflection_points_2trig = inflection_points_2trig_1hr
                    color_1trig = "blue"
                    color_2trig = "red"
                if len(inflection_points_1trig) >= 1:
                    mean_1trig = np.mean(inflection_points_1trig)
                    inflection_1 = plt.axvline(mean_1trig, color=color_1trig, linestyle="--")
                    inflection_handles_and_labels[inflection_1] = format_inflection_label("1Trig Inflection",mean_1trig)
                    legend_groups["group1_inflection"].append(inflection_1)
                if len(inflection_points_2trig) >= 1:
                    mean_2trig = np.mean(inflection_points_2trig)
                    inflection_2 = plt.axvline(mean_2trig, color=color_2trig, linestyle="--")
                    inflection_handles_and_labels[inflection_2] = format_inflection_label("2Trig Inflection",mean_2trig)
                    legend_groups["group2_inflection"].append(inflection_2)
                knockout_inflection_points = (
                    inflection_points_knockout_48hr
                    if any(test in tests for test in ["A280_48hr", "A280_72hr"])
                    else inflection_points_knockout_1hr
                )
                if knockout_inflection_points:
                    mean_knockout = np.mean(knockout_inflection_points)
                    knockout_inflection = plt.axvline(
                        mean_knockout, color="forestgreen", linestyle="--")
                    inflection_label = "Knockout Inflection"
                    if test_name in ["A280_1hr", "A280_48hr", "A280_72hr"]:
                        inflection_label = f"Knockout {time_label_from_test(test_name)} Inflection"
                    inflection_handles_and_labels[knockout_inflection] = format_inflection_label(
                        inflection_label, mean_knockout)
                    legend_groups["other"].append(knockout_inflection)
            if "A400" == tests[0]:
                y_axis = 'Normalized Absorbance Units'
                assay = 'Turbidity'
            else:
                y_axis = "Normalized Sedimentation"
                assay = 'Sedimentation'
            if len(protein_constructs) == 1:
                protein_title = format_chart_title(protein_constructs[0],len(protein_constructs))
                plt.title(f"{protein_title}, {assay}", fontsize = 14)
                plt.ylabel(y_axis)
                plt.xlabel("pH")
            else:
                if (len(protein_constructs) == 2 or is_alone_knockout_comparison) and len(tests) == 1:
                    # Use the first construct as the base protein title
                    protein_title = format_chart_title(protein_constructs[0],len(protein_constructs))
                    plt.title(f"{protein_title}, {assay}", fontsize = 14)
                    plt.ylabel(y_axis)
                    plt.xlabel("pH")
            # add asterisk legend entry LAST so it appears on bottom row
            if combined_any_omitted_points:
                    asterisk_handle = plt.scatter([], [], marker='$*$', color='black', s=40)
                    handles_and_labels[asterisk_handle] = "Omitted from curve fitting"
                    legend_groups["other"].append(asterisk_handle)
            handles_and_labels.update(inflection_handles_and_labels)
            ordered_handles = (legend_groups["group1_main"] + legend_groups["group1_inflection"] + legend_groups["group1_ci"] + legend_groups["group2_main"] + legend_groups["group2_inflection"] + legend_groups["group2_ci"] + legend_groups["other"])
            # remove duplicates while preserving order
            seen = set()
            handles = []
            for h in ordered_handles:
                if h not in seen and h in handles_and_labels:
                    handles.append(h)
                    seen.add(h)
            labels = [handles_and_labels[h] for h in handles]
            apply_custom_y_limits(protein_constructs, tests, plot_singly=False)
            ax = plt.gca()
            ax.yaxis.set_major_formatter(mticker.FormatStrFormatter('%.2f'))
            plt.legend(handles, labels,handler_map={tuple: HandlerTuple(ndivide=2, pad=1)},loc="upper right")
            safe_constructs = "_vs_".join([p.replace("/", "_") for p in protein_constructs])
            safe_tests = "_vs_".join(tests)
            plt.savefig(f"{BASE_DIR}/figures_with_CIs/{safe_constructs}__{safe_tests}.png", dpi=300, bbox_inches="tight")
            plt.show()
            plt.close()

if __name__ == '__main__':
    for protein_construct in protein_dict:
        if '2Trig' not in protein_construct and 'Gravity' not in protein_construct and 'Knockout' not in protein_construct:
            # Single vs. double-trigger A400
            if protein_construct == "1TEL Alone":
                statisticize([
                    ("1TEL Alone", "A400"),
                    ("2Trig-1TEL Alone", "A400"),
                    ("1TEL Knockout", "A400"),
                ], plot=True)
            else:
                statisticize([(f'{protein_construct}','A400'),(f'2Trig-{protein_construct}','A400')],plot=True)
    for protein_construct in protein_dict:
        if '2Trig' not in protein_construct and 'Gravity' not in protein_construct and 'Knockout' not in protein_construct:
            # Single vs. Double-trigger A280-1hr
            if protein_construct == "1TEL-PA-TNK1.UBA":
                statisticize([(f'{protein_construct}','A280_1hr'),(f'2Trig-{protein_construct}','A280_1hr')],plot=True,tailored=True)
            elif protein_construct == "1TEL Alone":
                statisticize([
                    ("1TEL Alone", "A280_1hr"),
                    ("2Trig-1TEL Alone", "A280_1hr"),
                    ("1TEL Knockout", "A280_1hr"),
                ], plot=True)
            else:
                statisticize([(f'{protein_construct}','A280_1hr'),(f'2Trig-{protein_construct}','A280_1hr')],plot=True)

    for protein_construct in protein_dict:
        if '2Trig' not in protein_construct and 'Gravity' not in protein_construct and 'Knockout' not in protein_construct:
            # Single-trigger A280 1hr vs 48-72 hr
            if protein_construct != '1TEL-PA-TNK1.UBA':
                statisticize([(f'{protein_construct}','A280_1hr'),(f'{protein_construct}','A280_48-72hr')],plot=True,tailored=True)
            else:
                statisticize([('1TEL-PA-TNK1.UBA', 'A280_1hr'),('1TEL-PA-TNK1.UBA', 'A280_48-72hr')],plot=True,tailored=True,specific_stock_id='S3')

        elif '2Trig' in protein_construct:
            # Double-trigger A280 1hr vs 48-72 hr
            statisticize([(f'{protein_construct}','A280_1hr'),(f'{protein_construct}','A280_48-72hr')],plot=True,tailored=True)
        else:
            # Single-trigger Gravity A280 1hr vs 48-72 hr
            statisticize([(f'{protein_construct}','A280_1hr'),(f'{protein_construct}','A280_48-72hr')],plot=True,tailored=True)

    # Single trigger TV-vWA A400, A280 1 HR, A280 48-72 HR, A280 1 HR vs 48-72 HR
    statisticize([('10xHis-1TEL-TV-vWA','A280_48-72hr'),('10xHis-1TEL-TV-vWA (Gravity)','A280_48-72hr')],plot=True)
    statisticize([('10xHis-1TEL-TV-vWA','A400'),('10xHis-1TEL-TV-vWA (Gravity)','A400')],plot=True)
    statisticize([('10xHis-1TEL-TV-vWA','A280_1hr'),('10xHis-1TEL-TV-vWA (Gravity)','A280_1hr')],plot=True)
    # Single-Trigger vs Double-Trigger 1TEL-PA-TNK1.UBA 1 Hour vs. 48 Hour with both behaviors
    statisticize([(f'1TEL-PA-TNK1.UBA','A280_1hr'), (f'2Trig-1TEL-PA-TNK1.UBA','A280_1hr')],plot=True,tailored=True)