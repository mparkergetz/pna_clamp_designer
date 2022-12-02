#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  5 15:31:23 2022

@author: getz
"""

import os
from datetime import datetime

name = '799f-1115f_region_110bp' # name of sequence to be blocked without file extension
db = 'current_Bacteria_unaligned'# name of database without file extension

is_pna = 'y' # is PNA? (y/n)
forward_primer_position = 'left' # 'left or 'right'

lower_limit = 15 # set blocker length limits
upper_limit = 17 # set blocker length limits

skip_blast = 'y'


########################### Generate blocker fasta ############################

input_file = open(f'../input/sequences/{name}.fasta','r')
input_seq_list = input_file.read().splitlines()
input_seq = input_seq_list[1]
seq_len = len(input_seq)
input_file.close()

in_path = '/home/busbyp/Desktop/blocker_gen/input'

## define PNA check functions ##

# function checks for strings of 6 or more purines
def purine_test(blocker_temp):
    pass_test = 0
    b = (len(blocker_temp))-6
    for x in range(0,b):
        count = 0
        for i in range(x,x+6):
            if blocker_temp[i] == "A" or blocker_temp[i] == "G":
                count +=1
            else:
                count = 0
        if count >= 6:
            pass_test += 1
    return (test + pass_test)

# function checks for < 50% purine content
def p50_test(blocker_temp):
    pass_test = 0
    p_count = 0
    for i in range(0,len(blocker_temp)):
        if blocker_temp[i] == "A" or blocker_temp[i] == "G":
            p_count += 1
    if (p_count/len(blocker_temp))>.55:
        pass_test = 1
    return (test + pass_test)

# function checks for <35% G content
def gcon_test(blocker_temp):    
    pass_test = 0
    g_count = 0
    for i in range(0,len(blocker_temp)):
        if blocker_temp[i] == "G":
            g_count += 1
    if (g_count/len(blocker_temp))>.38:
        pass_test = 1
    return (test + pass_test)

if is_pna == 'y':
    prefix = 'pna'
else: 
    prefix = 'blocker'
    
#set up output folder
out_path = f'/home/busbyp/Desktop/blocker_gen/output/{name}_{prefix}s'
try:
    os.mkdir(out_path)
    print('Output folder created')
except:
    print('Output folder exists')

#generate PNA dictionary, checking for parameters
blocker_dict = {}

test = 0 
pass_test = 0


name_no = 0
for l in range(lower_limit,upper_limit+1):
    x=0
    for i in range(0,(len(input_seq)-l)):
        blocker_name = f"{prefix}_{name_no+1}"
        blocker_temp = input_seq[x:(x+l)]
        if is_pna == 'y':
            if (purine_test(blocker_temp) + p50_test(blocker_temp) + gcon_test(blocker_temp) == 0):
                blocker_dict[blocker_name] = blocker_temp
                name_no +=1
        else:
            blocker_dict[blocker_name] = blocker_temp
            name_no +=1
        x += 1

# write out .fasta file
newline = '\n'
query_fasta = f'{out_path}/{prefix}s_for_{name}.fasta'
temp_fasta = open(query_fasta, "w")
for i in blocker_dict.keys():
    temp_fasta.write(f">{i}{newline}{blocker_dict[i]}{newline}")
temp_fasta.close()

# write out full list of names
keys = open(f"{out_path}/full_{prefix}_list_{name}.txt", "w+")
for i in blocker_dict.keys():
    keys.write(f"{i}{newline}")
#keys.close()

print('Blockers generated')


#### Set up blast db, run blastn based on specified parameters ###############

if os.path.isfile(f'{in_path}/blastdbs/'+db+'.nhr'):
    print('Database exists')
else:
    if os.path.isfile(f'{in_path}/sequences/'+db+'.fa'):
        print('Unzipped file present')
    else:
        os.system(f'gzip -dk {in_path}/sequences/{db}.fa.gz')
        print('Creating database')
        cmd = (f'makeblastdb -in {in_path}/sequences/{db}.fa '
              f'-dbtype nucl -title {db} '
              f'-out {in_path}/blastdbs/{db} '
              f'-logfile {in_path}/blastdbs/{db}.txt'
              )
        os.system(cmd)
    os.remove(f'{in_path}/sequences/{db}.fa')    
       
blast_args = ('-dust no ' 
              '-evalue 1000 ' 
              '-perc_identity 50 ' 
              '-word_size 8 '
              '-max_target_seqs 10000 '
              '-num_threads 12 '
              '-outfmt "10 qseqid sseqid pident qlen length evalue"'
              )

blocker_args = (f'-query {query_fasta} '
               f'-db {in_path}/blastdbs/{db} '
               f'-out {out_path}/{name}_{db}_output.csv'
               )


cmd = "blastn "+blast_args+" "+blocker_args
if skip_blast == 'n':
    now = datetime.now()
    current_time = now.strftime("%H:%M")
    print('Blastn start at',current_time)
    os.system(cmd) 
    now = datetime.now()
    current_time = now.strftime("%H:%M")
    print('Blastn end at',current_time)
else:
    print('Blast skipped')



###################### Filter results ########################################

import pandas as pd
from scipy import stats

# write in list of oligo clamp names
#keys = open("full_oligo_list.txt","r")
keys.seek(0)
oligo_out = keys.read()
oligo_out = oligo_out.rstrip('\n')
oligo_out_list = oligo_out.split("\n")

# self blast to determine positions
cmd = f"blastn -query {out_path}/{prefix}s_for_{name}.fasta -subject ../input/sequences/{name}.fasta -word_size {lower_limit} -evalue 1000 -dust no -perc_identity 100 -max_target_seqs 1 -out {out_path}/self_blast_output.csv -outfmt '10 qseqid sstart send'"
os.system(cmd)
colnames2 = ['qseqid', 'sstart', 'send']
len_df = pd.read_csv(f"{out_path}/self_blast_output.csv",names = colnames2, header = None )
if forward_primer_position == 'left':
    len_df.rename(columns = {'sstart':'dist'}, inplace = True)
else:
    len_df = len_df.assign(dist = (seq_len-len_df['send']))
os.remove(f"{out_path}/self_blast_output.csv")


# write in blastn output
colnames = ['qseqid', 'sseqid', 'pident', 'qlen', 'length', 'evalue']
df = pd.read_csv(out_path+f'/{name}_{db}_output.csv',names = colnames, header = None)

df_temp = df.assign(offset = (df['qlen']-df['length']))

'''
df_filter = df_temp.groupby('qseqid',as_index = False).min()
df_filter = df_filter[((df_filter['pident']>90) & (df_filter['offset'] < 1))]
df_pass = df_temp[~df_temp['qseqid'].isin(df_filter['qseqid'].tolist())]
'''
df_p100 = df_temp[((df_temp['pident']==100) & (df_temp['offset'] == 0))]
df_p100 = df_p100['qseqid'].value_counts()

df_counts = df_temp['qseqid'].value_counts()
df_counts = df_counts.reset_index()
df_counts.columns = ['qseqid','count']

df_temp = pd.merge(df_temp,len_df[['qseqid','dist']],on='qseqid', how='left')
df_temp = pd.merge(df_temp,df_counts[['qseqid','count']],on='qseqid', how='left')

df_scaled = df_temp
df_scaled[['pident','offset', 'dist', 'count']] = stats.zscore(df_scaled[['pident','offset', 'dist', 'count']])


df_scaled = df_scaled.groupby('qseqid')['pident','offset','dist','count'].agg('mean')
df_scaled = df_scaled.assign(score = ((-(df_scaled['offset']))+df_scaled['dist']+df_scaled['count']+df_scaled['pident']))
df_scaled = df_scaled.sort_values('score', ascending = True).reset_index(drop=False)

#df_final = df_scaled[['qseqid','score']]
#df_final = pd.merge(df_final,len_df[['qseqid','dist']],on='qseqid', how='left')

#df_final = df_final.loc[df_final.groupby('qseqid').score.idxmin()].reset_index(drop=True)
#df_final = df_final.sort_values('score', ascending = True).reset_index(drop=True)

#df_final.to_csv("./top_oligo_hits_standardized.csv")
