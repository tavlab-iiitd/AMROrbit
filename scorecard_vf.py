import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
import numpy as np
import os
import re
import random
import json

from adjustText import adjust_text
import csv
random.seed(10)

organisms = ['Escherichia coli', 'Klebsiella pneumoniae', 'Pseudomonas aeruginosa', 'Acinetobacter baumannii', 'Staphylococcus aureus', 'Enterococcus faecium', 'Enterobacter species']
sources = ['Blood', 'Urine', 'Wound', 'Sputum']


# organisms = ['Escherichia coli', 'Klebsiella pneumoniae', ]
# sources = ['Blood', 'Urine', 'Wound', 'Sputum']


for organism in organisms:
    for source in sources:
        print("Processing organism---------------------------------", organism)
        print("Processing source---------------------------------", source)
        # folders = ["ecoli_coli", "ecoli_imi", "ecoli_mero", "kleb_coli", "kleb_imi", "kleb_mero"]
        folders = os.listdir(f'Time series data/{organism}/{source}/')

        if '.DS_Store' in folders:
            folders.remove('.DS_Store')
        
        for folder in folders:
            antibiotic = folder.split('_')[0]
            path = os.path.join('Paper_Tables', organism, source)
            os.makedirs(path, exist_ok=True)
            csv_file = open(os.path.join(path, f'{organism}_{antibiotic}_{source}.csv'), "w")
            csv_writer = csv.writer(csv_file)
            row = ['Country', '2014-2017', '2015-2018', '2016-2019', '2017-2020', '2018-2021', '2019-2022']
            csv_writer.writerow(row)

            quadrant_data = {}

            print("Processing folder---------------------------------", folder)
            # Define folder path containing the files
            folder_path = 'Time series data/' + organism + '/' + source + '/' + folder

            # Regular expression to extract organism, antibiotic, and year from filenames
            file_pattern = re.compile(r"(?P<organism>[^_]+)_(?P<antibiotic>[^_]+)_I_(?P<year>\d{4})")

            # Create a list to store file paths and their corresponding metadata
            file_data = []

            # List all files in the folder
            for file_name in os.listdir(folder_path):
                match = file_pattern.match(file_name)
                if match:
                    metadata = match.groupdict()
                    metadata['file_path'] = os.path.join(folder_path, file_name)
                    file_data.append(metadata)
            if len(file_data) == 0:
                print(f"No files found in {folder_path}")
                continue
            # Group the files by organism and antibiotic
            grouped_files = {}
            for data in file_data:
                key = (data['organism'], data['antibiotic'])
                if key not in grouped_files:
                    grouped_files[key] = []
                grouped_files[key].append(data)

            # Country colors dictionary
            country_colors = {
                "Argentina": "blue", "Australia": "green", "Belgium": "lavender", "Brazil": "pink",
                "Canada": "orange", "Chile": "brown", "China": "purple", "Colombia": "grey",
                "Croatia": "cyan", "Czech Republic": "magenta", "Denmark": "yellow", "France": "lime",
                "Germany": "black", "Greece": "royalblue", "Hungary": "navy", "Ireland": "teal",
                "Israel": "olive", "Italy": "maroon", "Japan": "gold", "Korea, South": "lightblue",
                "Kuwait": "orchid", "Mexico": "peru", "Netherlands": "plum", "Panama": "salmon",
                "Philippines": "sienna", "Poland": "silver", "Portugal": "wheat", "South Africa": "tomato",
                "Spain": "turquoise", "Switzerland": "violet", "Taiwan": "tan", "Thailand": "indigo",
                "Turkey": "red", "United Kingdom": "crimson", "United States": "khaki", "Venezuela": "coral"
            }

            print("Grouped files---------------------------------", grouped_files)

            # Process each organism-antibiotic combination
            for (organism, antibiotic), files in grouped_files.items():
                # Determine global min and max for the current organism-antibiotic combination
                global_intercept_min = float('inf')
                global_intercept_max = float('-inf')
                global_slope_min = float('inf')
                global_slope_max = float('-inf')

                # Prepare to save correlation results
                correlation_results = []
                json_results = []
                for data in files:
                    df = pd.read_csv(data['file_path'])
                    assert(len(df) == len(df.dropna()))

                    country_data = df[df['Country'] != 'Global']
                    global_average = df[df['Country'] == 'Global'].iloc[0]

                    # Calculate new slope and intercept
                    country_data['slope'] += global_average['slope']
                    country_data['intercept'] += global_average['intercept']

                    # Shift intercepts to be positive and log transform
                    #country_data['intercept'] = np.log(country_data['intercept'] - country_data['intercept'].min() + 1)
                    #country_data['slope'] = np.log(country_data['slope'] - country_data['slope'].min() + 1)

                    # Update global min and max values
                    global_intercept_min = min(global_intercept_min, country_data['intercept'].min())
                    global_intercept_max = max(global_intercept_max, country_data['intercept'].max())
                    global_slope_min = min(global_slope_min, country_data['slope'].min())
                    global_slope_max = max(global_slope_max, country_data['slope'].max())

                    global_slope_max += 0.1
                    global_intercept_max += 0.1

                for data in files:
                    df = pd.read_csv(data['file_path'])

                    year = data['year']
                    correlation_results.append({data['year']})

                    assert(len(df) == len(df.dropna()))

                    country_data = df[df['Country'] != 'Global']
                    global_average = df[df['Country'] == 'Global'].iloc[0]

                    # Calculate new slope and intercept
                    country_data['slope'] += global_average['slope']
                    country_data['intercept'] += global_average['intercept']

                    #print("Before log-----------------------")
                    #print(country_data)

                    # Shift intercepts to be positive and log transform
                    #country_data['intercept'] = np.log(country_data['intercept'] - country_data['intercept'].min() + 1)
                    #country_data['slope'] = np.log(country_data['slope'] - country_data['slope'].min() + 1)


                    # Calculate median values for slope and intercept
                    median_slope = country_data['slope'].median()
                    median_intercept = country_data['intercept'].median()


                    # Adding jitter to avoid overlapping points
                    #jitter_strength = 0.02
                    #country_data['intercept'] += np.random.uniform(0, jitter_strength, size=len(country_data))
                    #country_data['slope'] += np.random.uniform(0, jitter_strength, size=len(country_data))

                    #print("After log-----------------------")
                    correlation_results.append(country_data)


                    plt.figure(figsize=(10, 8))
                    texts = []
                    for index, row in country_data.iterrows():
                        x = row['intercept']
                        y = row['slope']

                        if x < median_intercept and y < median_slope:
                            quadrant = 1

                        elif x >= median_intercept and y < median_slope:
                            quadrant = 2

                        elif x < median_intercept and y >= median_slope:
                            quadrant = 3

                        else:
                            quadrant = 4

                        country_name = row['Country']
                        if country_name not in quadrant_data:
                            quadrant_data[country_name] = {2014: 0, 2015: 0, 2016: 0, 2017: 0, 2018: 0, 2019: 0}
                        quadrant_data[country_name][int(data['year'])] = quadrant

                        if country_name in ['Taiwan']:
                            color = country_colors.get(country_name, 'black')
                            plt.scatter(row['intercept'], row['slope'], color=color, s=1000, label=country_name)

                        else:
                            color = country_colors.get(country_name, 'black')
                            plt.scatter(row['intercept'], row['slope'], color=color, s=120, label=country_name, alpha=0.3, linewidth=0.5)
                        #texts.append(plt.text(row['intercept'] + 0.03, row['slope'] - 0.05, country_name, fontsize=5, rotation=45))
                        # plt.text(row['intercept'] + 0.02, row['slope'] - 0.02, country_name, fontsize=4, rotation=30)

                    #adjust_text(texts, arrowprops=dict(arrowstyle="->", color='r', lw=0.5))

                    # Plot lines for median values with transparency
                    plt.axhline(y=median_slope, color='red', linestyle='--', label='Median Slope', alpha=0.5)
                    plt.axvline(x=median_intercept, color='green', linestyle='--', label='Median Intercept', alpha=0.5)

                    # Set global x and y limits
                    plt.xlim([global_intercept_min, global_intercept_max])
                    plt.ylim([global_slope_min, global_slope_max])


                    plt.title(f"{data['year']} - {int(data['year']) + 3}", fontsize=16)
                    plt.xlabel('Intercept', fontsize=12)
                    plt.ylabel('Slope', fontsize=12)
                    #plt.legend(loc='best', fontsize='xx-small', ncol=7)

                    # Ensure the directory for saving the figures exists
                    organism_dir = os.path.join('Plots', organism, source)
                    os.makedirs(organism_dir, exist_ok=True)

                    # Save the figure
                    plt.savefig(os.path.join(organism_dir, f"{antibiotic}_{data['year']}.pdf"), bbox_inches='tight')
                    plt.savefig(os.path.join(organism_dir, f"{antibiotic}_{data['year']}.png"), bbox_inches='tight')
                    plt.close()

                    plots_dir = os.path.join('Scorecards_plots', organism, source, antibiotic)
                    os.makedirs(plots_dir, exist_ok=True)

                    # Save the figure
                    plt.savefig(os.path.join(plots_dir, f"{antibiotic}_{data['year']}.png"), bbox_inches='tight')

                    # Correlation between slope and intercept
                    slope = country_data['slope']
                    intercept = country_data['intercept']


                    if len(slope) < 2 or len(intercept) < 2:
                        print(f"Not enough data points for correlation calculation for {organism} with {antibiotic} in {data['year']}. Skipping...")
                        continue
                    else:
                        correlation, p_value = pearsonr(slope, intercept)
                        print(f"Correlation between slope and intercept for {organism} with {antibiotic} in {data['year']}: {correlation:.4f} (p-value: {p_value:.4f})")

                    # Save the correlation result
                    correlation_results.append(f"Correlation between slope and intercept for {organism} with {antibiotic} in {data['year']}: {correlation:.4f} (p-value: {p_value:.4f})\n")

                    # Define the structure of the JSON
                    year_json= {
                        "year": year,
                        "median_slope": median_slope,
                        "median_intercept": median_intercept,
                        "countries": [
                            {
                                "name": row['Country'],
                                "x": row['intercept'],
                                "y": row['slope']
                            } for _, row in country_data.iterrows()
                        ]
                    }

                    json_results.append(year_json)
                    # Print JSON
                    print(json.dumps(year_json))

                f = open(os.path.join(organism_dir, (f"{antibiotic}_" + organism + "_year.json")), "w")
                # Convert to JSON
                json.dump(json_results, f, indent=4)
                f.close()

                # Save correlation results to a text file
                with open(os.path.join(organism_dir, f'{antibiotic}_correlation_results.txt'), 'w') as f:
                    f.writelines(str(correlation_results))

                for country in quadrant_data:
                    csv_writer.writerow([country, quadrant_data[country][2014], quadrant_data[country][2015], quadrant_data[country][2016], quadrant_data[country][2017], quadrant_data[country][2018], quadrant_data[country][2019]])

                

