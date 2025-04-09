import os
import pandas as pd
import csv


folder = 'Time series data'
infections = os.listdir(folder)
res = 0
count = 0

file = open('Median_shift_first.csv', 'w')
writer = csv.writer(file)

writer.writerow(['Infection', 'Source', 'Antibiotic', 'First subset intercept', 'Last subset intercpt', 'Shift'])

if '.DS_Store' in infections:
    infections.remove('.DS_Store')

for infection in infections:
    sources = os.listdir(f'Time series data/{infection}/')
    if '.DS_Store' in sources:
        sources.remove('.DS_Store')

    for source in sources:
        antibiotics = os.listdir(f'Time series data/{infection}/{source}')
        if '.DS_Store' in antibiotics:
            antibiotics.remove('.DS_Store')
        for antibiotic in antibiotics:
            files = os.listdir(f'Time series data/{infection}/{source}/{antibiotic}')

            if '.DS_Store' in files:
                files.remove('.DS_Store')

            files.sort()

            if not files:
                print(infection, source, antibiotic)
                continue
            
            file1_data = pd.read_csv(f'Time series data/{infection}/{source}/{antibiotic}/{files[0]}')
            file2_data = pd.read_csv(f'Time series data/{infection}/{source}/{antibiotic}/{files[-1]}')

            file1_intercept = file1_data.loc[file1_data['Country'] == 'Global', 'intercept'].iloc[0]
            file2_intercept =  file2_data.loc[file2_data['Country'] == 'Global', 'intercept'].iloc[0]

            

            if file2_intercept > file1_intercept:
                writer.writerow([infection, source, antibiotic, file1_intercept, file2_intercept, file2_intercept - file1_intercept])
                res+=1

            count+=1

print(res)
print(count)

res = 0
count = 0

file = open('Median_shift_continuous.csv', 'w')
writer = csv.writer(file)

writer.writerow(['Infection', 'Source', 'Antibiotic', '2014-2017', '2015-2018', '2016-2019', '2017-2020', '2018-2021', '2019-2022'])
if '.DS_Store' in infections:
    infections.remove('.DS_Store')

for infection in infections:
    sources = os.listdir(f'Time series data/{infection}/')
    if '.DS_Store' in sources:
        sources.remove('.DS_Store')

    for source in sources:
        antibiotics = os.listdir(f'Time series data/{infection}/{source}')
        if '.DS_Store' in antibiotics:
            antibiotics.remove('.DS_Store')
        for antibiotic in antibiotics:
            files = os.listdir(f'Time series data/{infection}/{source}/{antibiotic}')

            if '.DS_Store' in files:
                files.remove('.DS_Store')

            files.sort()

            if not files:
                print(infection, source, antibiotic)
                continue
            
            intercept = float('-inf')
            temp = {'2014-2017':None, '2015-2018':None, '2016-2019':None, '2017-2020':None, '2018-2021':None, '2019-2022':None}

            for file in files:

                data = pd.read_csv(f'Time series data/{infection}/{source}/{antibiotic}/{file}')
                intercept_ = data.loc[data['Country'] == 'Global', 'intercept'].iloc[0]

                if '2014_2017' in file:
                    temp['2014-2017'] = intercept_

                elif '2015_2018' in file:
                    temp['2015-2018'] = intercept_

                elif '2016_2019' in file:
                    temp['2016-2019'] = intercept_

                elif '2017_2020' in file:
                    temp['2017-2020'] = intercept_

                elif '2018_2021' in file:
                    temp['2018-2021'] = intercept_

                elif '2019_2022' in file:
                    temp['2019-2022'] = intercept_

                
                if intercept_ > intercept:
                    intercept = intercept_
                else:
                    break

            else:
                writer.writerow([infection, source, antibiotic] + list(temp.values()))
                res+=1

            count+=1

print(res)
print(count)