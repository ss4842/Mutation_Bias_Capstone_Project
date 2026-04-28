
import os
os.chdir("/Users/ss4842/Desktop")

#grabbing genomes
sample_names = []  #storing genome names
with open("ST_clean_singletons_test.vcf", "r") as f: #opening up my clean singletons vcf (recombination-free)
    while True: 
        line = f.readline()   #read vcf line by line
        if not line:
            break 
        
        if line.startswith("#CHROM"):     #pick up header rows
            sample_names = line.split()[9:]    #sample names = split columns 9 onwards in header row
            break 
            
sub_types = ["A>T", "A>C", "A>G", "T>A", "T>G", "T>C", "G>C", "G>A", "G>T", "C>G", "C>A", "C>T"] #all possible substitution types

sub_counts = {} #new dict to keep substitution counts
for i in range(len(sample_names)): #for each sample 
    sub_counts[i] = {} #nested dict for each sample to organise
    for sub in sub_types: #for each substitution type in list
        sub_counts[i][sub] = 0 #every sub type starts at 0 
      
with open("ST_clean_singletons_test.vcf", "r") as f: #opening vcf from the top
    while True:
        line = f.readline()   #read vcf line by line again
        if not line:
            break  
        
        if line.startswith("#") or len(line.strip()) == 0:  #skipping header lines & empty lines 
            continue
        
        entries = line.split("\t")
            
        if len(entries) < 10:  #if row has less than 10 columns, skips
            continue 
                    
        filter_status = entries[6]  #filtering for pass snps only
        if filter_status != "PASS":   
            continue 
                                       
        ref = entries[3]  #setting up the ref & alt labels to match all possible subtypes
        alt = entries[4] 
        sub_label = ref + ">" + alt 
                    
        if sub_label not in sub_types:  #skipping any label which isnt in sub_types list
            continue 
                                       
        for i, sample in enumerate(entries[9:]): #for each sample from column 9 onwards
            if sample == "1": #if there is a mutation 
                sub_counts[i][sub_label] += 1 #add 1 to the sample count for the particular sub label 

total_sub_counts = {}  #adding up all sub counts across all samples in new dict
for sub in sub_types:   #for each sub type in sub_types
    total_sub_counts[sub] = 0 #set count to 0 to start off
                                       
for i in range(len(sample_names)): #for each sample 
    for sub in sub_types:   #and for each sub type
        total_sub_counts[sub] += sub_counts[i][sub] #adding up total sub counts for each sample 
                                
grand_total = 0  #calculating total of all subs
for sub in sub_types:  #for each sub type
     grand_total += total_sub_counts[sub] #adding up all sub counts 
                                       
print("Overall substitution spectrum:")
print("-" * 40)
for sub in sub_types:
    print(sub, total_sub_counts[sub])

print("-" * 40)
print("Total:", grand_total)
                                       
with open("ST_substitution_counts_Final.csv", "w") as f:  #saving to csv
    header = "sample," + ",".join(sub_types)  #joining list of sub types into string separated by commas 
    f.write(header + "\n")  
    for i in range(len(sample_names)):  
        row = sample_names[i]    #writing genome names in rows
        for sub in sub_types:   
            row = row + "," + str(sub_counts[i][sub])  #adding sub count for each sub type in row
        f.write(row + "\n")


#substitution matrix
bases = ["A", "T", "G", "C"]   #the four valid bases
matrix = {}    #empty dict for matrix
for ref in bases:      #for each ref base
    matrix[ref] = {}    #inner dict created 
    for alt in bases:   #for each alt base
        if ref == alt:    #if ref and alt are same base, skip - no mut present
            continue
        matrix[ref][alt] = total_sub_counts[f"{ref}>{alt}"]    #counting for each ref>alt combo
        
with open("ST_substitution_count_matrix_Final.csv", "w") as f:
    f.write("REF/ALT,A,T,G,C\n")     #writing out header row
    for ref in bases:    #for each ref base
        row = ref      #row = ref base
        for alt in bases:    #for each alt base
            if ref == alt:   #if ref & alt = same base
                row = row + ",0"    #write out 0 (diagonally
            else:
                row = row + "," + str(matrix[ref][alt])   #every other ref>alt combo, add up count 
        f.write(row + "\n")   #write each row to csv
print("Matrix saved to ST_sub_count_matrix_V3.csv")

print("\nResults saved to substitution_counts_V3.csv")


#opportunity-normalisation
base_nos = {"C":0, "G":0, "A":0, "T":0} #starting all base counts at 0
with open("styphimurium.fna", "r") as f:   #opening reference fasta file
    for line in f:
        line = line.strip()     #removing white space
        if line.startswith(">"):     #skipping header line
            continue
        for base in line:         #going through each base in sequence
            if base in base_nos:     #if base in seq = valid base (skipping any bases that are not A, T, C, G)
                base_nos[base] += 1   #add 1 to base count for that particular base


for base in ["C", "G", "A", "T"]:     #printing no of each base in reference genome
    print(base, base_nos[base])   

#calculating denominators C+G & A+T
CG_denom = base_nos["C"] + base_nos["G"]   #calculating no of C bases + no of G bases to divide C-G subs by
AT_denom = base_nos["A"] + base_nos ["T"]  #calculating no of A bases + no of T bases to divide A-T subs by

print(CG_denom, AT_denom)      #checking CG content & AT content in reference genome

#6-class substitution spectra
CT = total_sub_counts["C>T"] + total_sub_counts["G>A"]   #C>T/G>A same substitution on diff strands, so add together
CG = total_sub_counts["C>G"] + total_sub_counts["G>C"]   #C>G/G>C same sub on diff strands - add together
CA = total_sub_counts["C>A"] + total_sub_counts["G>T"]   #C>A/G>T same sub on diff strands - add together 
TA = total_sub_counts["T>A"] + total_sub_counts["A>T"]   #T>A/A>T same sub on diff strands - add together
TC = total_sub_counts["T>C"] + total_sub_counts["A>G"]   #T>C/A>G same sub on diff strands - add together
TG = total_sub_counts["T>G"] + total_sub_counts["A>C"]   #T>G/A>C same sub on diff strands - add together

normalised_CT = CT/CG_denom    #opp-normalised - total C>T+G>A counts divided by total CG content 
normalised_CG = CG/CG_denom    #opp-normalised - total C>G+G>C counts divided by total CG content 
normalised_CA = CA/CG_denom    #opp-normalised - total C>A+G>T counts divided by total CG content 
normalised_TA = TA/AT_denom    #opp-normalised - total T>A+A>T counts divided by total AT content 
normalised_TC = TC/AT_denom    #opp-normalised - total T>C+A>G counts divided by total AT content 
normalised_TG = TG/AT_denom    #opp-normalised - total T>G+A>C counts divided by total AT content

#re-normalising for comparison on graph
total = normalised_CT + normalised_CG + normalised_CA + normalised_TA + normalised_TC + normalised_TG  #total of all normalised values

#dividing each normalised value by the total to scale values (so they add up to one)
renormalised_CT = normalised_CT/total
renormalised_CG = normalised_CG/total
renormalised_CA = normalised_CA/total
renormalised_TA = normalised_TA/total
renormalised_TC = normalised_TC/total
renormalised_TG = normalised_TG/total

with open("ST_Normalised_sub_table_final.csv", "w") as f:   #writing out in csv
    f.write("substitution,normalised_count,per_million_bases, renormalised_count\n")  #writing out headers
##for each of 6-class sub types, writing out the normalised value, what that would look like per 1mil bases, then the renormalised value to base graph off of
##each normalised value written to 6sf, per mil bases to 2sf, renormalised to 4sf
    f.write(f"C>T/G>A, {normalised_CT:.6f}, {normalised_CT * 1000000:.2f}, {renormalised_CT:.4f}\n")  
    f.write(f"C>G/G>C, {normalised_CG:.6f}, {normalised_CG * 1000000:.2f}, {renormalised_CG:.4f}\n")
    f.write(f"C>A/G>T, {normalised_CA:.6f}, {normalised_CA * 1000000:.2f}, {renormalised_CA:.4f}\n")
    f.write(f"T>A/A>T, {normalised_TA:.6f}, {normalised_TA * 1000000:.2f}, {renormalised_TA:.4f}\n")
    f.write(f"T>C/A>G, {normalised_TC:.6f}, {normalised_TC * 1000000:.2f}, {renormalised_TC:.4f}\n")
    f.write(f"T>G/A>C, {normalised_TG:.6f}, {normalised_TG * 1000000:.2f}, {renormalised_TG:.4f}\n")

print("Normalised substitutions saved to ST_normalised_sub_table_V3.csv")

#counting out transitions & transversions + calculating Ts/Tv ratio 
transitions = CT + TC    #transitions = C>T/G>A + T>C/A>G
transversions = CA + CG + TA + TG    #transversions = C>A/G>T + C>G/G>C + T>A/A>T + T>G/A>C
ts_tv = transitions/transversions   #Ts/Tv ratio = total transition count divided by total transversion count

print(f"Transitions : {transitions}")
print(f"Transversions: {transversions}")
print(f"Ts/Tv ratio: {ts_tv:.2f}")

with open("ST_Ts_Tv_ratio.csv", "w") as f:   #saving so csv
    f.write("species, transitions, transversions, ts_tv_ratio\n")    #writing out headers
    f.write(f"S. typhimurium, {transitions}, {transversions}, {ts_tv:.2f}\n")    #writing out columns

print("Ts/Tv ratio saved to ST_Ts_Tv_ratio.csv")




