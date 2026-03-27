import csv
from scipy.optimize import curve_fit
from scipy.stats import t
import matplotlib.pyplot as plt
from matplotlib.legend_handler import HandlerTuple
import numpy as np
from pathlib import Path

#1)Pull specified constructs and their corresponding specified data from csv and put into protein_dict by calling the master function 'statisticize'. Ex.:
    #statisticize([('Construct1','A400'),('Construct2','A400')])
#2)Repeat steps 3-4 for each construct asked for.
#3)For the purposes of data analysis, only combine biological replicates right before making plots.
#4)Determine medians, try to fit a curve, and plot.

#First, place all the data within a dictionary, separating by protein construct, then by date/concentration/pH, then by test-type.
protein_dict = {}

BASE_DIR = Path(__file__).resolve().parent
CSV_PATH = BASE_DIR / "data.csv"

with open(CSV_PATH, newline="", encoding="utf-8-sig") as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        protein_construct = row['Protein Construct'].strip()
        if protein_construct not in protein_dict:
            protein_dict[protein_construct] = {}
        bio_rep_key = row['Date'] + "_" + row['Stock Concentration mg/mL'] + '_mg/mL_pH_' + row['pH']
        if bio_rep_key not in protein_dict[protein_construct]:
            protein_dict[protein_construct][bio_rep_key] = {'A400': [],'A280_1hr': [],'A280_48-72hr': [],'Stock ID': row.get('Stock ID', '').strip(),'vWA ID': row.get('vWA ID', '').strip()}
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
        if row['A280_48-72hr']!="":
            raw_value = row['A280_48-72hr'].strip()
            if '*' in raw_value or '^' in raw_value or '~' in raw_value:
                protein_dict[protein_construct][row['Date']+"_"+row['Stock Concentration mg/mL']+'_mg/mL_pH_'+row['pH']]['A280_48-72hr'].append(raw_value)
            else:
                protein_dict[protein_construct][row['Date']+"_"+row['Stock Concentration mg/mL']+'_mg/mL_pH_'+row['pH']]['A280_48-72hr'].append(float(raw_value))

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
    is_solubility = test_name in ["A280_1hr", "A280_48-72hr"]
    # Case 1: 1Trig vs 2Trig
    if len(protein_constructs) == 2 and len(tests) == 1:
        if is_vwa and any("Gravity" in c for c in protein_constructs):
            base_label = "Gravity" if is_gravity else "ÄKTA"
        else:
            base_label = "2Trig" if is_2trig else "1Trig"
        if test_name == "A280_1hr":
            label = f"{base_label} 1 Hour"
        elif test_name == "A280_48-72hr":
            label = f"{base_label} 48 Hour"
        else:
            label = base_label
    # Case 2: 1 hr vs 48 hr for same construct
    elif len(protein_constructs) == 1 and has_1hr and has_48hr:
        if is_gravity:
            trig_label = "Gravity"
        else:
            trig_label = "2Trig" if is_2trig else "1Trig"
        time_label = "1 Hour" if test_name == "A280_1hr" else "48 Hour"
        label = f"{trig_label} {time_label}"
    # Case 3: single plot, one construct and one test
    elif len(protein_constructs) == 1 and len(tests) == 1:
        if is_gravity:
            base_label = "Gravity"
        else:
            base_label = "2Trig" if is_2trig else "1Trig"
        if is_solubility:
            time_label = "1 Hour" if test_name == "A280_1hr" else "48 Hour"
            label = f"{base_label} {time_label}"
        else:
            label = base_label
    else:
        label = construct
    return f"{label} (n = {n})"

def display_inflection_x(inflection_value):
    if inflection_value <= 4.0:
        return 4.5
    elif inflection_value >= 9.0:
        return 8.5
    else:
        return inflection_value

def format_inflection_label(prefix, inflection_value):
    if inflection_value <= 4.0:
        return f"{prefix}: pH ≤ 4.5"
    elif inflection_value >= 9.0:
        return f"{prefix}: pH ≥ 8.5"
    else:
        return f"{prefix}: pH = {inflection_value:.1f}"

def get_dates_for_construct(protein_construct, test):
    if protein_construct not in protein_dict:
        return []
    dates = set()
    for bio_rep in protein_dict[protein_construct]:
        date = bio_rep.split("_")[0]
        values = protein_dict[protein_construct][bio_rep][test]
        if len(values) != 0:
            dates.add(date)
    return sorted(dates)

def behavior_label_from_vwa_id(vwa_id):
    if vwa_id == "W1":
        return "Behavior 1"
    elif vwa_id == "W2":
        return "Behavior 2"
    return vwa_id

def get_vwa_behavior_colors(base_color):
    return {"W1": base_color,"W2": tuple(max(0, c * 0.65) for c in base_color[:3]) + (base_color[3],) if len(base_color) == 4
              else tuple(max(0, c * 0.65) for c in base_color[:3])}

def apply_custom_y_limits(protein_constructs, tests, plot_singly, specific_date=None):
    key = (tuple(sorted(protein_constructs)), tuple(sorted(tests)))
    single_plot_limits = {}
    combined_plot_limits = {
        (tuple(sorted(("1TEL-PA-TNK1.UBA", "2Trig-1TEL-PA-TNK1.UBA"))), ("A400",)): {"ymax": 5.0},
        (tuple(sorted(("10xHis-1TEL-TV-vWA", "10xHis-1TEL-TV-vWA (Gravity)"))), ("A280_1hr",)): {"ymin": -1.0},
        (("1TEL-GG-TNK1.UBA",), ("A280_1hr", "A280_48-72hr")): {"ymin": -1.0},
        (("2Trig-1TEL-GG-TNK1.UBA",), ("A280_1hr", "A280_48-72hr")): {"ymin": -0.2},
        (("1TEL-PA-TNK1.UBA",), ("A280_1hr", "A280_48-72hr")): {"ymin": -0.4},
        (tuple(sorted(("2Trig-10xHis-1TEL-TV-vWA", "10xHis-1TEL-TV-vWA"))), ("A280_1hr",)): {"ymin": -0.8},
        (tuple(sorted(("2Trig-1TEL-GG-TNK1.UBA", "1TEL-GG-TNK1.UBA"))), ("A280_1hr",)): {"ymin": -0.4},
        (tuple(sorted(("2Trig-1TEL-PA-TNK1.UBA", "1TEL-PA-TNK1.UBA"))), ("A280_1hr",)): {"ymin": -0.5},
        (tuple(sorted(("2Trig-10xHis-1TEL-SR-TNK1.UBA", "10xHis-1TEL-SR-TNK1.UBA"))), ("A400",)): {"ymin": -0.2, "ymax": 1.8},}
    limits = single_plot_limits.get(key) if plot_singly else combined_plot_limits.get(key)
    if limits is not None:
        plt.ylim(bottom=limits.get("ymin"), top=limits.get("ymax"))

def statisticize(proteins_and_tests:list,plot_singly=False,plot=False,tailored=False,pool=False, specific_date=None,specific_stock_id=None):
    """
    Normally, you should pass only one test and one protein construct into 'proteins_and_tests', in the form of a list of tuples, where the
    tuple looks like: (protein_construct,test)
    'protein_constructs' are the protein_dict dictionary keys that refer to the particular protein constructs of interest.
    'tests' are the types of tests you performed, referred to by the variable being measured, or the column name in the csv. Ie.: 'A400', 'A280_1hr'
    """
    #The first half of this function is just a filter pulling applicable data from the protein_dict dictionary.
    data = {'protein construct':[],'date':[],'pH':[],'test':[],'data':[],'tailored_data':[], 'is_starred':[], 'is_single_omitted':[], 'is_time_omitted':[], 'vwa_id':[]}
    protein_constructs = []
    tests = []
    for protein_construct,test in proteins_and_tests:
        if protein_construct not in protein_constructs:
            protein_constructs.append(protein_construct)
        if test not in tests:
            tests.append(test)
    # Detect 1 hr vs 48 hr comparison
    is_time_comparison = (len(protein_constructs) == 1 and "A280_1hr" in tests and "A280_48-72hr" in tests)
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
                has_48hr = len(protein_dict[protein_construct][bio_rep]["A280_48-72hr"]) != 0
                if not (has_1hr and has_48hr):
                    continue
            if len(protein_dict[protein_construct][bio_rep][test])!=0:
                #if the biological replicate actually had this test done, aggregate its data into the dataframe replicate by replicate, conserving their order.
                for replicate in protein_dict[protein_construct][bio_rep][test]:
                    data['protein construct'].append(protein_construct)
                    data['date'].append(date)
                    pH = bio_rep.split("_")[-1]
                    data['pH'].append(pH)
                    data['test'].append(test)
                    vwa_id = protein_dict[protein_construct][bio_rep].get('vWA ID', '').strip()
                    data['vwa_id'].append(vwa_id)

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
                    if test in ['A280_1hr', 'A280_48-72hr']:
                        value = (theoretical_max_concentration - numeric_replicate) / theoretical_max_concentration
                    # Normalized Absorbance: Abs @ 400 nm * (Expected Concentration/Actual Concentration)
                    else:
                        value = numeric_replicate * (5/concentration)
                    data['data'].append(value)
                    data['tailored_data'].append(value)
                    data['is_starred'].append(is_starred)
                    data['is_single_omitted'].append(is_single_omitted)
                    data['is_time_omitted'].append(is_time_omitted)

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
        is_single_omitted = data['is_single_omitted'][replicate]
        is_time_omitted = data['is_time_omitted'][replicate]
        if pH not in constructs_by_pH_by_date[construct_key]:
            constructs_by_pH_by_date[construct_key][pH] = {"dates": {},"dates_starred": {},"dates_single_omitted": {},"dates_time_omitted": {},"dates_vwa_id": {}}
        if date not in constructs_by_pH_by_date[construct_key][pH]["dates"]:
            constructs_by_pH_by_date[construct_key][pH]["dates"][date] = []
            constructs_by_pH_by_date[construct_key][pH]["dates_starred"][date] = []
            constructs_by_pH_by_date[construct_key][pH]["dates_single_omitted"][date] = []
            constructs_by_pH_by_date[construct_key][pH]["dates_time_omitted"][date] = []
        constructs_by_pH_by_date[construct_key][pH]["dates"][date].append(measurement)
        constructs_by_pH_by_date[construct_key][pH]["dates_starred"][date].append(is_starred)
        constructs_by_pH_by_date[construct_key][pH]["dates_single_omitted"][date].append(is_single_omitted)
        constructs_by_pH_by_date[construct_key][pH]["dates_time_omitted"][date].append(is_time_omitted)
        constructs_by_pH_by_date[construct_key][pH]["dates_vwa_id"][date] = data['vwa_id'][replicate]

    biological_n_by_construct = {}

    for construct in constructs_by_pH_by_date:
        unique_dates = set()
        for pH in constructs_by_pH_by_date[construct]:
            for date in constructs_by_pH_by_date[construct][pH]["dates"]:
                unique_dates.add(date)
        biological_n_by_construct[construct] = len(unique_dates)

    biological_n_by_construct_and_vwa = {}

    for construct in constructs_by_pH_by_date:
        biological_n_by_construct_and_vwa[construct] = {}
        for pH in constructs_by_pH_by_date[construct]:
            for date in constructs_by_pH_by_date[construct][pH]["dates"]:
                vwa_id = constructs_by_pH_by_date[construct][pH]["dates_vwa_id"].get(date, "").strip()
                if vwa_id not in biological_n_by_construct_and_vwa[construct]:
                    biological_n_by_construct_and_vwa[construct][vwa_id] = set()
                biological_n_by_construct_and_vwa[construct][vwa_id].add(date)
    # convert sets to counts
    for construct in biological_n_by_construct_and_vwa:
        for vwa_id in biological_n_by_construct_and_vwa[construct]:
            biological_n_by_construct_and_vwa[construct][vwa_id] = len(biological_n_by_construct_and_vwa[construct][vwa_id])

    # Compute the 95% CI Band for each pH value within each construct
    alpha_bonferroni = 0.05/11
    for construct in constructs_by_pH_by_date:
        for pH in constructs_by_pH_by_date[construct]:
            pH_dates = constructs_by_pH_by_date[construct][pH]["dates"]
            pH_starred = constructs_by_pH_by_date[construct][pH]["dates_starred"]
            pH_single_omitted = constructs_by_pH_by_date[construct][pH]["dates_single_omitted"]
            pH_time_omitted = constructs_by_pH_by_date[construct][pH]["dates_time_omitted"]

            all_drops_for_display = []
            all_drops_for_fit = []

            for date_key in pH_dates:
                drops = pH_dates[date_key]
                star_flags = pH_starred[date_key]
                single_omit_flags = pH_single_omitted[date_key]
                time_omit_flags = pH_time_omitted[date_key]

                paired = list(zip(drops, star_flags, single_omit_flags, time_omit_flags))

                # Data shown in CI band:
                # include starred points, but still honor single-plot-only omissions
                display_paired = paired

                # Data used for fitting:
                fit_paired = paired
                if tailored:
                    fit_paired = [(d, s, so, to) for d, s, so, to in fit_paired if not s]
                if plot_singly:
                    fit_paired = [(d, s, so, to) for d, s, so, to in fit_paired if not so]
                if is_time_comparison:
                    fit_paired = [(d, s, so, to) for d, s, so, to in fit_paired if not to]

                display_drops = [d for d, s, so, to in display_paired]
                fit_drops = [d for d, s, so, to in fit_paired]

                if len(display_drops) > 0:
                    all_drops_for_display.extend(np.array(display_drops, dtype=float))
                if len(fit_drops) > 0:
                    all_drops_for_fit.extend(np.array(fit_drops, dtype=float))

            # display stats for CI
            if len(all_drops_for_display) > 0:
                N_display = len(all_drops_for_display)
                sd_display = np.std(all_drops_for_display)
                SE_display = sd_display / np.sqrt(N_display)
                t_display = t.ppf(1 - (alpha_bonferroni / 2), df=N_display - 1)
                error_display = SE_display * t_display
                mean_display = np.mean(all_drops_for_display)

                constructs_by_pH_by_date[construct][pH]["mean_display"] = mean_display
                constructs_by_pH_by_date[construct][pH]["CI_lower_display"] = mean_display - error_display
                constructs_by_pH_by_date[construct][pH]["CI_upper_display"] = mean_display + error_display

            # fit stats for curve fitting
            if len(all_drops_for_fit) > 0:
                N_fit = len(all_drops_for_fit)
                sd_fit = np.std(all_drops_for_fit)
                SE_fit = sd_fit / np.sqrt(N_fit)
                t_fit = t.ppf(1 - (alpha_bonferroni / 2), df=N_fit - 1)
                error_fit = SE_fit * t_fit
                mean_fit = np.mean(all_drops_for_fit)

                constructs_by_pH_by_date[construct][pH]["mean"] = mean_fit
                constructs_by_pH_by_date[construct][pH]["CI_lower"] = mean_fit - error_fit
                constructs_by_pH_by_date[construct][pH]["CI_upper"] = mean_fit + error_fit
                constructs_by_pH_by_date[construct][pH]["N"] = N_fit

    def pH_to_absorbance_model_4pl(pH,upper_asymptote,Hill_slope,inflection_point,lower_asymptote):
        #This is the model (4 Parameter Logistic, 4PL) that we're going to try to fit using scipy's curve_fit function.
        return lower_asymptote + (upper_asymptote - lower_asymptote) / (1 + (pH / inflection_point)**Hill_slope)

    def classify_legend_group(construct, test_name, protein_constructs, tests):
        is_gravity_comparison = (len(protein_constructs) == 2 and len(tests) == 1 and any("Gravity" in c for c in protein_constructs))
        is_time_comparison = (len(protein_constructs) == 1 and "A280_1hr" in tests and "A280_48-72hr" in tests)

        if is_gravity_comparison:
            if "Gravity" in construct:
                return "group2"
            return "group1"  # ÄKTA first

        if is_time_comparison:
            if test_name == "A280_48-72hr":
                return "group2"
            return "group1"  # 1 Hour first

        # default: 1Trig vs 2Trig
        if "2Trig" in construct:
            return "group2"
        return "group1"

    if plot==True:
        plt.figure(figsize=(8.75,6))
        handles_and_labels = {}
        inflection_handles_and_labels = {}
        legend_groups = {"group1_main": [],"group1_inflection": [],"group1_ci": [],"group2_main": [],"group2_inflection": [],"group2_ci": [],"other": []}
        is_combined_vwa_1hr_solubility = (len(protein_constructs) == 2 and len(tests) == 1 and tests[0] == "A280_1hr" and "10xHis-1TEL-TV-vWA" in protein_constructs)
        combined_any_omitted_points = False

        def change_colors():
            keys_1trig = []
            keys_2trig = []
            keys_48_hour = []
            keys_2Trig_48_hour = []
            keys_Gravity = []
            keys_Gravity_48_hour = []
            for construct_key in constructs_by_pH_by_date:
                construct, test_name = construct_key
                if "Gravity" in construct and test_name == "A280_48-72hr":
                    keys_Gravity_48_hour.append(construct_key)
                elif "Gravity" in construct:
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
            brown_shades = shade_list(plt.cm.copper, len(keys_Gravity_48_hour))
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
            for construct, color in zip(sorted(keys_Gravity_48_hour), brown_shades):
                color_map[construct] = color
            return color_map
        color_map = change_colors()

        inflection_points_1trig_1hr = []
        inflection_points_1trig_48hr = []
        inflection_points_2trig_1hr = []
        inflection_points_2trig_48hr = []
        inflection_points_gravity_1hr = []
        inflection_points_gravity_48hr = []

        for construct_key in constructs_by_pH_by_date:
            any_starred_points = False
            construct, test_name = construct_key
            is_special_vwa_construct = (is_combined_vwa_1hr_solubility and test_name == "A280_1hr" and "vWA" in construct and "2Trig" not in construct and "Gravity" not in construct)
            pH_values_fit = []
            means_fit = []
            pH_values_display = []
            ci_lower_bounds = []
            ci_upper_bounds = []
            construct_pHs =[]
            construct_data = []
            construct_starred = []
            construct_single_omitted = []
            construct_vwa_ids = []
            construct_time_omitted = []
            special_curve_handles = {}
            special_scatter_handles = {}
            # sort pHs numerically
            sorted_pHs = sorted(constructs_by_pH_by_date[construct_key].keys(), key=float)
            vwa_groups = {}
            for pH in sorted_pHs:
                pH_entry = constructs_by_pH_by_date[construct_key][pH]
                if is_special_vwa_construct:
                    for date in pH_entry["dates"]:
                        vwa_id = pH_entry["dates_vwa_id"].get(date, "").strip()
                        if vwa_id not in vwa_groups:
                            vwa_groups[vwa_id] = {'pH_values': [], 'means': [], 'ci_lower': [], 'ci_upper': []}
                # Always keep raw points for plotting
                for date, measurements in pH_entry["dates"].items():
                    star_flags = pH_entry["dates_starred"][date]
                    single_omit_flags = pH_entry["dates_single_omitted"][date]
                    time_omit_flags = pH_entry["dates_time_omitted"][date]
                    vwa_id = pH_entry["dates_vwa_id"].get(date, "").strip()
                    for measurement, is_starred, is_single_omitted, is_time_omitted in zip(measurements, star_flags, single_omit_flags, time_omit_flags):
                        construct_pHs.append(float(pH))
                        construct_data.append(measurement)
                        construct_starred.append(is_starred)
                        construct_single_omitted.append(is_single_omitted)
                        construct_vwa_ids.append(vwa_id)
                        construct_time_omitted.append(is_time_omitted)
                # Only use pHs with computed statistics for fitting / CI bands
                # Skip only if there are no stats at all for either display or fitting
                if "mean" not in pH_entry and "mean_display" not in pH_entry:
                    continue
                if is_special_vwa_construct:
                    subgroup_measurements = {}
                    for date, measurements in pH_entry["dates"].items():
                        vwa_id = pH_entry["dates_vwa_id"].get(date, "").strip()
                        if vwa_id not in subgroup_measurements:
                            subgroup_measurements[vwa_id] = []
                        star_flags = pH_entry["dates_starred"][date]
                        single_omit_flags = pH_entry["dates_single_omitted"][date]
                        time_omit_flags = pH_entry["dates_time_omitted"][date]
                        paired = list(zip(measurements, star_flags, single_omit_flags, time_omit_flags))
                        if tailored:
                            paired = [(d, s, so, to) for d, s, so, to in paired if not s]
                        if plot_singly:
                            paired = [(d, s, so, to) for d, s, so, to in paired if not so]
                        if is_time_comparison:
                            paired = [(d, s, so, to) for d, s, so, to in paired if not to]
                        subgroup_measurements[vwa_id].extend([d for d, s, so, to in paired])
                    for vwa_id, vals in subgroup_measurements.items():
                        if len(vals) == 0:
                            continue
                        vals = np.array(vals, dtype=float)
                        mean_val = np.mean(vals)
                        if len(vals) > 1:
                            sd_val = np.std(vals, ddof=1)
                            se_val = sd_val / np.sqrt(len(vals))
                            t_val = t.ppf(1 - (alpha_bonferroni / 2), df=len(vals) - 1)
                            err = se_val * t_val
                        else:
                            err = 0.0
                        vwa_groups[vwa_id]['pH_values'].append(float(pH))
                        vwa_groups[vwa_id]['means'].append(mean_val)
                        vwa_groups[vwa_id]['ci_lower'].append(mean_val - err)
                        vwa_groups[vwa_id]['ci_upper'].append(mean_val + err)
                else:
                    # Add CI-band values whenever display stats exist
                    if "mean_display" in pH_entry:
                        pH_values_display.append(float(pH))
                        ci_lower_bounds.append(pH_entry["CI_lower_display"])
                        ci_upper_bounds.append(pH_entry["CI_upper_display"])

                    # Add fit values only when fit stats exist
                    if "mean" in pH_entry:
                        pH_values_fit.append(float(pH))
                        means_fit.append(pH_entry["mean"])

            try:
                if is_special_vwa_construct:
                    for vwa_id, group_data in vwa_groups.items():
                        if len(group_data['pH_values']) < 3:
                            continue
                        group_pH = group_data['pH_values']
                        group_means = group_data['means']
                        group_ci_lower = group_data['ci_lower']
                        group_ci_upper = group_data['ci_upper']
                        pH_linspace = np.linspace(min(group_pH), max(group_pH), 400)
                        initial_guesses = [max(group_means), -4, 5.5, min(group_means)]
                        best_fit_parameters, covariance_matrix = curve_fit(pH_to_absorbance_model_4pl,group_pH,group_means,p0=initial_guesses)
                        inflection_value = best_fit_parameters[2]
                        if inflection_value <= 4.0:
                            continue
                        vwa_colors = get_vwa_behavior_colors(color_map[construct_key])
                        curve_handle_list = plt.plot(pH_linspace,pH_to_absorbance_model_4pl(pH_linspace, *best_fit_parameters),color=vwa_colors[vwa_id],linestyle="-" if vwa_id == "W1" else "-")
                        line_handle = curve_handle_list[0]
                        ci_handle = plt.fill_between(group_pH,group_ci_lower,group_ci_upper,color=vwa_colors[vwa_id],alpha=0.12)
                        handles_and_labels[ci_handle] = f"{behavior_label_from_vwa_id(vwa_id)} 95% CI"
                        legend_groups["group1_ci"].append(ci_handle)
                        inflection_x = display_inflection_x(inflection_value)
                        inflection_handle = plt.axvline(inflection_x,color=vwa_colors[vwa_id],linestyle=":")
                        inflection_handles_and_labels[inflection_handle] = (
                            f"{behavior_label_from_vwa_id(vwa_id)} Inflection: "
                            f"{format_inflection_label('', inflection_value).replace(': ', '').strip()}")
                        legend_groups["group1_inflection"].append(inflection_handle)
                        special_curve_handles[vwa_id] = line_handle
                else:
                    pH_linspace = np.linspace(min(pH_values_fit), max(pH_values_fit), 400)
                    initial_guesses = [max(means_fit), -4, 5.5, min(means_fit)]
                    best_fit_parameters, covariance_matrix = curve_fit(pH_to_absorbance_model_4pl, pH_values_fit, means_fit, p0=initial_guesses)
                    inflection_value = best_fit_parameters[2]
                    if inflection_value <= 4.0:
                        raise ValueError(f"Calculated inflection point ({inflection_value:.1f}) is <= 4.0")
                    if "Gravity" in construct:
                        if test_name == "A280_48-72hr":
                            inflection_points_gravity_48hr.append(inflection_value)
                        else:
                            inflection_points_gravity_1hr.append(inflection_value)
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
                    curve_handle_list = plt.plot(pH_linspace,pH_to_absorbance_model_4pl(pH_linspace, *best_fit_parameters),color=color_map[construct_key])
                    line_handle = curve_handle_list[0]
                    if "Gravity" in construct and test_name == 'A280_48-72hr':
                        ci_handle = plt.fill_between(pH_values_display, ci_lower_bounds, ci_upper_bounds, color='brown', alpha=0.20)
                        ci_label = "Gravity 48 Hour 95% CI"
                    elif "Gravity" in construct and test_name == 'A280_1hr':
                        ci_handle = plt.fill_between(pH_values_display, ci_lower_bounds, ci_upper_bounds, color='lime', alpha=0.20)
                        ci_label = "Gravity 1 Hour 95% CI"
                    elif "Gravity" in construct:
                        ci_handle = plt.fill_between(pH_values_display, ci_lower_bounds, ci_upper_bounds, color='lime', alpha=0.20)
                        ci_label = "Gravity 95% CI"
                    elif test_name == "A280_48-72hr" and "2Trig" not in construct:
                        ci_handle = plt.fill_between(pH_values_display, ci_lower_bounds, ci_upper_bounds, color='pink', alpha=0.20)
                        ci_label = "1Trig 48 Hour 95% CI"
                    elif test_name == "A280_48-72hr" and "2Trig" in construct:
                        ci_handle = plt.fill_between(pH_values_display, ci_lower_bounds, ci_upper_bounds, color='mediumpurple', alpha=0.20)
                        ci_label = "2Trig 48 Hour 95% CI"
                    elif "2Trig" in construct:
                        ci_handle = plt.fill_between(pH_values_display, ci_lower_bounds, ci_upper_bounds, color='tomato', alpha=0.20)
                        ci_label = "2Trig 95% CI"
                    else:
                        ci_handle = plt.fill_between(pH_values_display, ci_lower_bounds, ci_upper_bounds, color='cyan', alpha=0.20)
                        ci_label = "1Trig 95% CI"
                    handles_and_labels[ci_handle] = ci_label
                    group_prefix = classify_legend_group(construct, test_name, protein_constructs, tests)
                    legend_groups[f"{group_prefix}_ci"].append(ci_handle)

                # plot raw scatter points
                jitter_strength = 0.08
                jittered_pHs = np.array(construct_pHs) + np.random.uniform(-jitter_strength, jitter_strength,len(construct_pHs))
                if is_special_vwa_construct:
                    special_scatter_handles = {}
                    vwa_colors = get_vwa_behavior_colors(color_map[construct_key])
                    for vwa_id, marker_style in [("W1", "^"), ("W2", "s")]:
                        xs = [x for x, gid in zip(jittered_pHs, construct_vwa_ids) if gid == vwa_id]
                        ys = [y for y, gid in zip(construct_data, construct_vwa_ids) if gid == vwa_id]
                        if len(xs) > 0:
                            scatter_handle = plt.scatter(xs,ys,marker=marker_style,color=vwa_colors[vwa_id],alpha=0.2)
                            special_scatter_handles[vwa_id] = scatter_handle
                    for vwa_id in special_curve_handles:
                        line_handle = special_curve_handles[vwa_id]
                        scatter_handle = special_scatter_handles.get(vwa_id)
                        behavior_label = behavior_label_from_vwa_id(vwa_id)
                        subgroup_n = biological_n_by_construct_and_vwa.get(construct_key, {}).get(vwa_id, 0)
                        base_label = format_legend_label(construct,test_name,protein_constructs,tests,subgroup_n)
                        if scatter_handle is not None:
                            pair_handle = (line_handle, scatter_handle)
                            handles_and_labels[pair_handle] = f"{base_label}, {behavior_label}"
                            legend_groups["group1_main"].append(pair_handle)
                        else:
                            handles_and_labels[line_handle] = f"{base_label}, {behavior_label}"
                            legend_groups["group1_main"].append(line_handle)
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
                if not is_special_vwa_construct:
                    legend_label = format_legend_label(construct,test_name,protein_constructs,tests,biological_n_by_construct[construct_key])
                    pair_handle = (line_handle, scatter_handle)
                    handles_and_labels[pair_handle] = legend_label
                    group_prefix = classify_legend_group(construct, test_name, protein_constructs, tests)
                    legend_groups[f"{group_prefix}_main"].append(pair_handle)

                ###--------------------------------------------------------------------------###

                if plot_singly:
                    inflection_x = display_inflection_x(inflection_value)
                    if "Gravity" in construct and test_name == "A280_48-72hr":
                        inflection_color = 'brown'
                    elif "Gravity" in construct:
                        inflection_color = "green"
                    elif test_name == "A280_48-72hr" and "2Trig" in construct:
                        inflection_color = "indigo"
                    elif test_name == "A280_48-72hr" and "2Trig" not in construct:
                        inflection_color = "purple"
                    elif "2Trig" in construct:
                        inflection_color = "red"
                    else:
                        inflection_color = "blue"
                    inflection_handle = plt.axvline(inflection_x,color=inflection_color,linestyle="--")
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
                if "Gravity" in construct:
                    ci_handle = plt.fill_between(pH_values_display, ci_lower_bounds, ci_upper_bounds, color='lime', alpha=0.20)
                elif test_name == "A280_48-72hr" and "2Trig" not in construct:
                    ci_handle = plt.fill_between(pH_values_display, ci_lower_bounds, ci_upper_bounds, color='pink', alpha=0.20)
                elif test_name == "A280_48-72hr" and "2Trig" in construct:
                    ci_handle = plt.fill_between(pH_values_display, ci_lower_bounds, ci_upper_bounds, color='mediumpurple', alpha=0.20)
                elif "2Trig" in construct:
                    ci_handle = plt.fill_between(pH_values_display, ci_lower_bounds, ci_upper_bounds, color='tomato', alpha=0.20)
                else:
                    ci_handle = plt.fill_between(pH_values_display, ci_lower_bounds, ci_upper_bounds, color='cyan', alpha=0.20)
                jitter_strength = 0.08
                jittered_pHs = np.array(construct_pHs) + np.random.uniform(-jitter_strength, jitter_strength, len(construct_pHs))
                if "2Trig" in construct and test_name == 'A280_48-72hr':
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
                handles_and_labels[ci_handle] = f"{legend_label.split(' (n')[0]} 95% CI"
                group_prefix = classify_legend_group(construct, test_name, protein_constructs, tests)
                legend_groups[f"{group_prefix}_ci"].append(ci_handle)
                legend_groups[f"{group_prefix}_main"].append(scatter_handle)
                if plot_singly and construct == "2Trig-10xHis-1TEL-SR-DARPin":
                    if test_name == "A280_48-72hr":
                        inflection_color = "indigo"
                    else:
                        inflection_color = "red"
                    inflection_handle = plt.axvline(4.5, color=inflection_color, linestyle="--")
                    handles_and_labels[inflection_handle] = "Inflection Point: pH ≤ 4.5"
                if (plot_singly and construct == "10xHis-1TEL-TV-vWA" and test_name == "A400" and specific_date == "1/28/2026"):
                    inflection_handle = plt.axvline(4.5, color="blue", linestyle="--")
                    handles_and_labels[inflection_handle] = "Inflection Point: pH ≤ 4.5"

            # mean points with Bonferroni error bars
            #plt.errorbar(x=pH_values,y=means, yerr=[np.array(means) - np.array(ci_lower_bounds), np.array(ci_upper_bounds) - np.array(means)], fmt='none', ecolor=color_map[construct], capsize=4, alpha=0.9)
            ###--------------------------------------------------------------------------###

            if plot_singly:
                if test_name == "A400":
                    y_axis = "Normalized Absorbance Units"
                    assay = "Turbidity"
                else:
                    y_axis = "Precipitation Fractional Loss"
                    assay = "Solubility"
                replicate_dates = sorted({d
                    for pH in constructs_by_pH_by_date[construct_key]
                    for d in constructs_by_pH_by_date[construct_key][pH]["dates"]})
                date_text = replicate_dates[0] if len(replicate_dates) == 1 else ", ".join(replicate_dates)
                plt.title(f"{format_protein_title(construct)}, {assay}\nReplicate Date: {date_text}",fontsize=14)
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
                apply_custom_y_limits(protein_constructs, tests, plot_singly=True, specific_date=specific_date)
                plt.legend(handles,labels,handler_map={tuple: HandlerTuple(ndivide=2, pad=1)},loc="upper right")
                safe_date_text = date_text.replace("/", "-").replace(", ", "_").replace(" ", "_")
                plt.savefig(f"{construct.replace('/', '_')}_{test_name}_{safe_date_text}_single.png",dpi=300,bbox_inches="tight")
                plt.show()
                plt.close()
                handles_and_labels = {}

        if not plot_singly:
            is_gravity_comparison = (len(protein_constructs) == 2 and len(tests) == 1 and any("Gravity" in construct for construct in protein_constructs))
            is_time_comparison = (len(protein_constructs) == 1 and "A280_1hr" in tests and "A280_48-72hr" in tests)
            if is_time_comparison:
                if "2Trig" in protein_constructs[0]:
                    if len(inflection_points_2trig_1hr) >= 1:
                        mean_1hr = np.mean(inflection_points_2trig_1hr)
                        inflection_x_1hr = display_inflection_x(mean_1hr)
                        inflection_1hr = plt.axvline(inflection_x_1hr, color="darkorange", linestyle="--")
                        inflection_handles_and_labels[inflection_1hr] = format_inflection_label("1 Hour Inflection",mean_1hr)
                        legend_groups["group1_inflection"].append(inflection_1hr)
                    elif protein_constructs[0] == "2Trig-10xHis-1TEL-SR-DARPin":
                        inflection_1hr = plt.axvline(4.5, color="darkorange", linestyle="--")
                        inflection_handles_and_labels[inflection_1hr] = "1 Hour Inflection: pH ≤ 4.5"
                        legend_groups["group1_inflection"].append(inflection_1hr)

                    if len(inflection_points_2trig_48hr) >= 1:
                        mean_48hr = np.mean(inflection_points_2trig_48hr)
                        inflection_x_48hr = display_inflection_x(mean_48hr)
                        inflection_48hr = plt.axvline(inflection_x_48hr, color="indigo", linestyle="--")
                        inflection_handles_and_labels[inflection_48hr] = format_inflection_label("48 Hour Inflection",mean_48hr)
                        legend_groups["group2_inflection"].append(inflection_48hr)
                    elif protein_constructs[0] == "2Trig-10xHis-1TEL-SR-DARPin":
                        inflection_48hr = plt.axvline(4.5, color="purple", linestyle="--")
                        inflection_handles_and_labels[inflection_48hr] = "48 Hour Inflection: pH ≤ 4.5"
                        legend_groups["group2_inflection"].append(inflection_48hr)

                elif "Gravity" in protein_constructs[0]:
                    if len(inflection_points_gravity_1hr) >= 1:
                        mean_1hr = np.mean(inflection_points_gravity_1hr)
                        inflection_x_1hr = display_inflection_x(mean_1hr)
                        inflection_1hr = plt.axvline(inflection_x_1hr, color="green", linestyle="--")
                        inflection_handles_and_labels[inflection_1hr] = format_inflection_label("1 Hour Inflection", mean_1hr)
                        legend_groups["group1_inflection"].append(inflection_1hr)

                    if len(inflection_points_gravity_48hr) >= 1:
                        mean_48hr = np.mean(inflection_points_gravity_48hr)
                        inflection_x_48hr = display_inflection_x(mean_48hr)
                        inflection_48hr = plt.axvline(inflection_x_48hr, color="brown", linestyle="--")
                        inflection_handles_and_labels[inflection_48hr] = format_inflection_label("48 Hour Inflection", mean_48hr)
                        legend_groups["group2_inflection"].append(inflection_48hr)   
                
                else:
                    if len(inflection_points_1trig_1hr) >= 1:
                        mean_1hr = np.mean(inflection_points_1trig_1hr)
                        inflection_x_1hr = display_inflection_x(mean_1hr)
                        inflection_1hr = plt.axvline(inflection_x_1hr, color="blue", linestyle="--")
                        inflection_handles_and_labels[inflection_1hr] = format_inflection_label("1 Hour Inflection",mean_1hr)
                        legend_groups["group1_inflection"].append(inflection_1hr)

                    if len(inflection_points_1trig_48hr) >= 1:
                        mean_48hr = np.mean(inflection_points_1trig_48hr)
                        inflection_x_48hr = display_inflection_x(mean_48hr)
                        inflection_48hr = plt.axvline(inflection_x_48hr, color="purple", linestyle="--")
                        inflection_handles_and_labels[inflection_48hr] = format_inflection_label("48 Hour Inflection",mean_48hr)
                        legend_groups["group2_inflection"].append(inflection_48hr)
            elif is_gravity_comparison:
                if "A400" in tests:
                    nongravity_points = inflection_points_1trig_1hr
                    nongravity_color = "blue"
                    nongravity_label = "ÄKTA Inflection"

                    gravity_points = inflection_points_gravity_1hr
                    gravity_color = "green"
                    gravity_label = "Gravity Inflection"

                elif "A280_48-72hr" in tests:
                    nongravity_points = inflection_points_1trig_48hr
                    nongravity_color = "purple"
                    nongravity_label = "ÄKTA 48 Hour Inflection"

                    gravity_points = inflection_points_gravity_48hr
                    gravity_color = "brown"
                    gravity_label = "Gravity 48 Hour Inflection"

                elif "A280_1hr" in tests:
                    nongravity_points = inflection_points_1trig_1hr
                    nongravity_color = "blue"
                    nongravity_label = "ÄKTA 1 Hour Inflection"

                    gravity_points = inflection_points_gravity_1hr
                    gravity_color = "green"
                    gravity_label = "Gravity 1 Hour Inflection"

                if len(nongravity_points) >= 1:
                    mean_ng = np.mean(nongravity_points)
                    inflection_x_ng = display_inflection_x(mean_ng)
                    inflection_ng = plt.axvline(inflection_x_ng, color=nongravity_color, linestyle="--")
                    inflection_handles_and_labels[inflection_ng] = format_inflection_label(nongravity_label, mean_ng)
                    legend_groups["group1_inflection"].append(inflection_ng)

                if len(gravity_points) >= 1:
                    mean_g = np.mean(gravity_points)
                    inflection_x_g = display_inflection_x(mean_g)
                    inflection_g = plt.axvline(inflection_x_g, color=gravity_color, linestyle="--")
                    inflection_handles_and_labels[inflection_g] = format_inflection_label(gravity_label, mean_g)
                    legend_groups["group2_inflection"].append(inflection_g)
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
                    inflection_x_1 = display_inflection_x(mean_1trig)
                    inflection_1 = plt.axvline(inflection_x_1, color=color_1trig, linestyle="--")
                    inflection_handles_and_labels[inflection_1] = format_inflection_label("1Trig Inflection",mean_1trig)
                    legend_groups["group1_inflection"].append(inflection_1)
                if len(inflection_points_2trig) >= 1:
                    mean_2trig = np.mean(inflection_points_2trig)
                    inflection_x_2 = display_inflection_x(mean_2trig)
                    inflection_2 = plt.axvline(inflection_x_2, color=color_2trig, linestyle="--")
                    inflection_handles_and_labels[inflection_2] = format_inflection_label("2Trig Inflection",mean_2trig)
                    legend_groups["group2_inflection"].append(inflection_2)
                else:
                    is_darpin_comparison = (
                            len(protein_constructs) == 2
                            and any(c == "10xHis-1TEL-SR-DARPin" for c in protein_constructs)
                            and any(c == "2Trig-10xHis-1TEL-SR-DARPin" for c in protein_constructs))
                    if is_darpin_comparison:
                        inflection_2 = plt.axvline(4.50, color=color_2trig, linestyle="--")
                        inflection_handles_and_labels[inflection_2] = "2Trig Inflection: pH ≤ 4.5"
                        legend_groups["group2_inflection"].append(inflection_2)
            if "A400" == tests[0]:
                y_axis = 'Normalized Absorbance Units'
                assay = 'Turbidity'
            else:
                y_axis = "Precipitation Fractional Loss"
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
            apply_custom_y_limits(protein_constructs, tests, plot_singly=False, specific_date=specific_date)
            plt.legend(handles, labels,handler_map={tuple: HandlerTuple(ndivide=2, pad=1)},loc="upper right")
            safe_constructs = "_vs_".join([p.replace("/", "_") for p in protein_constructs])
            safe_tests = "_vs_".join(tests)
            plt.savefig(f"{safe_constructs}__{safe_tests}.png", dpi=300, bbox_inches="tight")
            plt.show()
            plt.close()

if __name__ == '__main__':
    for protein_construct in protein_dict:
        if '2Trig' not in protein_construct and 'Gravity' not in protein_construct:
            # Single-trigger A400: one plot per replicate date
            for date in get_dates_for_construct(f'{protein_construct}', 'A400'):
                statisticize([(f'{protein_construct}', 'A400')],plot=True,plot_singly=True,specific_date=date)
            # Double-trigger A400: one plot per replicate date
            for date in get_dates_for_construct(f'2Trig-{protein_construct}', 'A400'):
                statisticize([(f'2Trig-{protein_construct}', 'A400')],plot=True,plot_singly=True,specific_date=date)
            for date in get_dates_for_construct('10xHis-1TEL-TV-vWA (Gravity)', 'A400'):
                statisticize([('10xHis-1TEL-TV-vWA (Gravity)', 'A400')],plot=True, plot_singly=True, specific_date=date)
            # Single vs. Double-trigger A400
            statisticize([(f'{protein_construct}','A400'),(f'2Trig-{protein_construct}','A400')],plot=True)

    for protein_construct in protein_dict:
        if '2Trig' not in protein_construct and 'Gravity' not in protein_construct:
            # Single-trigger A280 1 hr: one plot per replicate date
            for date in get_dates_for_construct(f'{protein_construct}', 'A280_1hr'):
                statisticize([(f'{protein_construct}', 'A280_1hr')],plot=True,plot_singly=True,specific_date=date)
            # Double-trigger A280 1 hr: one plot per replicate date
            for date in get_dates_for_construct(f'2Trig-{protein_construct}', 'A280_1hr'):
                statisticize([(f'2Trig-{protein_construct}', 'A280_1hr')],plot=True,plot_singly=True,specific_date=date)
            # Single vs. Double-trigger A280-1hr
            statisticize([(f'{protein_construct}','A280_1hr'),(f'2Trig-{protein_construct}','A280_1hr')],plot=True)

    for protein_construct in protein_dict:
        if '2Trig' not in protein_construct and 'Gravity' not in protein_construct:
            # Single-trigger A280 48-72 hr: one plot per replicate date
            for date in get_dates_for_construct(f'{protein_construct}', 'A280_48-72hr'):
               statisticize([(f'{protein_construct}', 'A280_48-72hr')],plot=True,plot_singly=True,specific_date=date)
            # Double-trigger A280 48-72 hr: one plot per replicate date
            for date in get_dates_for_construct(f'2Trig-{protein_construct}', 'A280_48-72hr'):
                statisticize([(f'2Trig-{protein_construct}', 'A280_48-72hr')],plot=True,plot_singly=True,specific_date=date)
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

    # Plot the individual Gravity replicates for A400, A280 1 HR, and A280 48-72 HR
    for date in get_dates_for_construct('10xHis-1TEL-TV-vWA (Gravity)', 'A400'):
        statisticize([('10xHis-1TEL-TV-vWA (Gravity)', 'A400')], plot=True, plot_singly=True, specific_date=date)
    for date in get_dates_for_construct('10xHis-1TEL-TV-vWA (Gravity)', 'A280_1hr'):
        statisticize([('10xHis-1TEL-TV-vWA (Gravity)', 'A280_1hr')], plot=True, plot_singly=True, specific_date=date)
    for date in get_dates_for_construct('10xHis-1TEL-TV-vWA (Gravity)', 'A280_48-72hr'):
        statisticize([('10xHis-1TEL-TV-vWA (Gravity)', 'A280_48-72hr')], plot=True, plot_singly=True, specific_date=date)
    # Single trigger TV-vWA A400, A280 1 HR, A280 48-72 HR, A280 1 HR vs 48-72 HR
    statisticize([('10xHis-1TEL-TV-vWA','A280_48-72hr'),('10xHis-1TEL-TV-vWA (Gravity)','A280_48-72hr')],plot=True)
    statisticize([('10xHis-1TEL-TV-vWA','A400'),('10xHis-1TEL-TV-vWA (Gravity)','A400')],plot=True)
    statisticize([('10xHis-1TEL-TV-vWA','A280_1hr'),('10xHis-1TEL-TV-vWA (Gravity)','A280_1hr')],plot=True)