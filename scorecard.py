import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
import numpy as np
import os
import re
import random
import json

from adjustText import adjust_text

random.seed(10)

folders = ["ecoli_coli", "ecoli_imi", "ecoli_mero", "kleb_coli", "kleb_imi", "kleb_mero"]


for folder in folders:
    print("Processing folder---------------------------------", folder)
    # Define folder path containing the files
    folder_path = '../processed_data/blood/' + folder + '/'

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

    # Group the files by organism and antibiotic
    grouped_files = {}
    for data in file_data:
        key = (data['organism'], data['antibiotic'])
        if key not in grouped_files:
            grouped_files[key] = []
        grouped_files[key].append(data)

    # Country colors dictionary
    country_colors = {
        "Argentina": "blue", "Australia": "green", "Belgium": "lavender", "Brazil": "purple",
        "Canada": "orange", "Chile": "brown", "China": "pink", "Colombia": "grey",
        "Croatia": "cyan", "Czech Republic": "magenta", "Denmark": "yellow", "France": "royalblue",
        "Germany": "black", "Greece": "lime", "Hungary": "navy", "Ireland": "teal",
        "Israel": "olive", "Italy": "maroon", "Japan": "gold", "Korea, South": "lightblue",
        "Kuwait": "orchid", "Mexico": "peru", "Netherlands": "plum", "Panama": "salmon",
        "Philippines": "sienna", "Poland": "silver", "Portugal": "tan", "South Africa": "tomato",
        "Spain": "turquoise", "Switzerland": "violet", "Taiwan": "wheat", "Thailand": "indigo",
        "Turkey": "red", "United Kingdom": "khaki", "United States": "crimson", "Venezuela": "coral"
    }

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

            # Calculate median values for slope and intercept
            median_slope = country_data['slope'].median()
            median_intercept = country_data['intercept'].median()

            print(country_data)


            plt.figure(figsize=(10, 8))
            texts = []
            for index, row in country_data.iterrows():
                country_name = row['Country']
                color = country_colors.get(country_name, 'black')
                plt.scatter(row['intercept'], row['slope'], color=color, s=300)#, label=country_name)
                #texts.append(plt.text(row['intercept'] + 0.03, row['slope'] - 0.05, country_name, fontsize=10, rotation=45))

                #plt.text(row['intercept'] + 0.02, row['slope'] - 0.02, country_name, fontsize=4, rotation=30)


            # Plot lines for median values with transparency
            plt.axhline(y=median_slope, color='red', linestyle='--', label='Median Slope', alpha=0.5, linewidth=5)
            plt.axvline(x=median_intercept, color='green', linestyle='--', label='Median Intercept', alpha=0.5, linewidth=5)

            # Set global x and y limits
            plt.xlim([global_intercept_min, global_intercept_max])
            plt.ylim([global_slope_min, global_slope_max])


            plt.title(f"{data['year']} - {antibiotic}")
            plt.xlabel('Intercept')
            plt.ylabel('Slope')
            plt.legend()

            # Ensure the directory for saving the figures exists
            organism_dir = os.path.join('../plots/', organism)
            os.makedirs(organism_dir, exist_ok=True)

            # Save the figure
            plt.savefig(os.path.join(organism_dir, f"{antibiotic}_{data['year']}.pdf"), bbox_inches='tight')
            plt.savefig(os.path.join(organism_dir, f"{antibiotic}_{data['year']}.png"), bbox_inches='tight')
            plt.close()


