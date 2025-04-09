import pandas
import os
import pandas as pd
import csv

organisms = ['Escherichia coli', 'Klebsiella pneumoniae', 'Pseudomonas aeruginosa', 'Acinetobacter baumannii', 'Staphylococcus aureus', 'Enterococcus faecium', 'Enterobacter species']
sources = ['Blood', 'Urine', 'Wound', 'Sputum']

csv_file = open('organism_antibiotic_source_trajectories.csv', 'w')
csv_writer = csv.writer(csv_file)
csv_writer.writerow(['Organism', 'Source', 'Antibiotic', 'Spiral In', 'Spiral Out', 'Constant', 'Total', 'Spiral In Count', 'Spiral Out Count', 'Constant Count', 'Spiral In %', 'Spiral Out %', 'Constant %'])



for organism in organisms:    
    for source in sources:
        files = os.listdir(f'Paper_Tables/{organism}/{source}/')

        for file in files:

            data = pd.read_csv(f'Paper_Tables/{organism}/{source}/{file}')
            spiral_in = []
            spiral_out = []
            constant = []
            total = 0

            for index, row in data.iterrows():
                total += 1
                country_name = row['Country']
                start = row['2014-2017']
                end = row['2019-2022']

                if start == end:
                    constant.append(country_name)

                elif end == 4 and start in [1,2,3]:
                    spiral_out.append(country_name)

                elif end == 1 and start in [2,3,4]:
                    spiral_in.append(country_name)

            print(file)
            print('Spiral In:')
            for country in spiral_in:
                print(country)

            print()

            print('Spiral Out:')
            for country in spiral_out:
                print(country)

            print()
            print('Constant:')
            for country in constant:
                print(country)
            
            print()
            print('Total Countries:', total)
            print('***********************************')

            l = file.split('_')
            antibiotic = l[1]

            if total != 0:
                csv_writer.writerow([organism, source, antibiotic, ','.join(spiral_in), ','.join(spiral_out), ','.join(constant), total, len(spiral_in), len(spiral_out), len(constant), len(spiral_in)/total * 100, len(spiral_out)/total * 100, len(constant)/total * 100])

csv_file.close()


data = pd.read_csv('organism_antibiotic_source_trajectories.csv')

new_csv = open('organism_antibiotic_source_count.csv', 'w')
csv_writer = csv.writer(new_csv)
csv_writer.writerow(['Organism', 'Source', 'Spiral In Count', 'Spiral Out Count', 'Constant Count', 'Total', 'Spiral In %', 'Spiral Out %', 'Constant %'])
for organism in organisms:
    for source in sources:
        rows = data[(data['Organism'] == organism) & (data['Source'] == source)]
        spiral_in_count = 0
        spiral_out_count = 0
        constant_count = 0
        total = 0
        for index, row in rows.iterrows():
            spiral_in_count += row['Spiral In Count']
            spiral_out_count += row['Spiral Out Count']
            constant_count += row['Constant Count']
            total += row['Total']
        
        print(organism, source, len(rows))
        csv_writer.writerow([organism, source, spiral_in_count, spiral_out_count, constant_count, total, spiral_in_count/total * 100, spiral_out_count/total * 100, constant_count/total * 100])

new_csv.close()