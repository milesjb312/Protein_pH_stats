'statisticize' can be considered the "main" function because this is where we do all of our calculations and plotting. Normally, you should pass only one test and one protein construct 
into 'proteins_and_tests', in the form of a list of tuples, where the tuple looks like: (protein_construct,test). 'protein_constructs' are the protein_dict dictionary keys 
that refer to the particular protein constructs of interest. 'tests' are the types of tests you performed, referred to by the variable being measured, or the column name in 
the csv (ie: 'A400', 'A280_1hr')

First, we make a new dictionary ("data") that we will pull most of our data from for calculations and plotting. We have several labels (eg 'is_starred', 'is_single_omitted',
'is_time_omitted','vwa_id') that we will use to identify data that needs to be omitted for curve fitting. All data points, including those labeled, are used to construct the 
95% confidence intervals. We use additional labels (eg "is_time_comparison") for graphing purposes. As part of our data reorganization we also calculate and store the theoretical
protein concentration that each eppendorf tube contained after adding protein to the pH buffer (this is our theoretical_max_concentration variable). We calculate it by dividing the
known protein stock concentration by five since we diluted the protein sample 5-fold (8 µL Stock Protein + 32 µL pH buffer). We normalize all of the solubility data by dividing the
difference between the theoretical max concentration and the actual measurement value by the theoretical max concentration. The resulting value is unitless and represents the fraction
of protein that precipitates out of solution, hence why we label the y-axis for solubility graphs as "Precipitation Fractional Loss" and why these are called solubility graphs since
we are reporting the percent of remaining soluble protein. To clarify, we use the theoretical max concentration tied to the replicate we are analyzing, which helps account for 
variability between protein stock concentrations. For the turbidity plots, we report the Normalized Absorbance Units by taking the A400 absorbance reading and  
dividing by 5 and the protein stock concentration to get the turbidity per mg/mL protein.

Below are brief explanations of the labels used in "statisticize":

is_time_comparison = Used to determine if the plot is comparing the 1 Hour data to it's 48 Hour counterpart. This should eventually be replaced with a logic check instead.
stock_id = Used to determine the 48 Hour counterpart to the 1 Hour data (there is a "Stock ID" column in the csv that we use to assign the stock IDs).
has_1hr = The replicate we're analyzing has 1 Hour solubility data. This should eventually be consolidated into a Stock ID logic check.
has_48hr = The replicate we're analyzing has 48 Hour solubility data. This should eventually be consolidated into a Stock ID logic check.
vwa_id = Used to determine which vWA behavior we observe in this set of data; this only applies to vWA solubility data.

Second, we iterate through each replicate and do the curve fitting and plotting. At the start, we need to define variables that will make it easier to do the curve fitting and plotting 
(eg "pH = data['pH'][replicate]). For calculating the 95% confidence intervals, we use the bonferroni method to minimize the chance of getting a type I error. We calculate a t-statistic
using the Student's t-distribution, although we don't believe there is enough statistical evidence to show that our data follow a Student's t-distribution. For the purposes of this study,
we are making the assumption that our data follow a Student's t-distribution purely for visualization purposes.

Third, we prepare the data for plotting and fit a four-parameter logistic (4PL) equation to which ever dataset we're interested in (eg all of the replicates for the 1Trig 1TEL-GG-TNK1.UBA turbidity assays). This uses parameters calculated using scipy's curve fitting function (command: curve_fit()). Scipy uses non-linear least squares to fit a function, f, to data. The actual curve fitting only takes up a few lines of code. Most of the code at this point deals with creating the actual plot and formating labels, colors, etc. The code tries fitting a 4PL curve to a dataset (single or combined as specified in the statisticize command). After fitting a curve, the code plots the fitted curve, 95% confidence interval as a shaded band, the raw data points, and the inflection-point line. A large portion of the code deals with making the legend labels. There is also a portion that initially fits two separate curves to the different vWA datasets that show unique behaviors. If curve fitting fails, then the code plots only the raw data points and in some instances includes an inflection-point line if we previously determined that it needed one (eg 2Trig-10xHis-1TEL-[sr]-DARPin). The code also ensures that the plot includes asterisks for any data points that were omitted from curve fitting. 
