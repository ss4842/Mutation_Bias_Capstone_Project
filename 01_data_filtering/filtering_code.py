import os
os.chdir("/Users/ss4842/Desktop/")

sample_names = []
with open("ST_raw.vcf", "r") as f: #opening vcf file
    while True:
        line = f.readline()  #reading line by line
        if not line:
            break
        if line.startswith("#CHROM"):    #focusing on header lines
            sample_names = line.split()[10:]  #splitting into columns - adding columns 10+ (skipping reference) to sample_names (the 1000 genomes)
            break       #once columns finish, break

sample_singleton_pos = {}  #nested dict 
for sample in sample_names: 
    sample_singleton_pos[sample] = {} #empty dict for each of the genomes - will store contig name + position of singletons

only_bases = {"A", "T", "G", "C"}  #only valid bases 

with open("ST_raw.vcf", "r") as f:  #opening raw vcf file
    for line in f:
        if line.startswith("#"):  #skipping headers
            continue
        entries = line.strip().split("\t")  #splitting into columns
        if len(entries) <10:  #skipping lines w less than 10 columns (would mean something is wrong there)
            continue
        if entries[6] != "PASS":     #filtering for "PASS"-only SNPs
            continue
        ref = entries[3]     #reference base = column 3(4)
        alt = entries[4]     #mutated base = column 4(5)
        if len(ref) != 1 or len(alt) != 1:  #skipping rows if they have more than one character i.e. multi-allelic sites
            continue
        if ref not in only_bases or alt not in only_bases:   #skipping rows if character present is not one of the four valid bases
            continue


#isolating singletons
        sample_genotypes = entries[10:]  #genotype values for every genome; value either "0" or "1"
        carriers = sum(1 for g in sample_genotypes if g == "1") #goes through every genotype in sample_genotypes - everytime a 1 is present, adds 1 to "carriers" for each genome
        if carriers != 1: #if number of genomes carrying particular mutation is not 1, meaning more than one genome has the mutation,
            continue         #not singleton so skip row
        contig = entries[0] #storing contigs
        pos = int(entries[1])  #turning pos into integer for recomb later
        for i, g in enumerate(sample_genotypes):   #looking at each genotype in list + their carrier value
            if g == "1":     #if carrier value = 1, meaning only 1 genome posseses particular mutation,
                sample = sample_names[i]      #looking for genome name in sample_list & storing in sample
                if contig not in sample_singleton_pos[sample]:   #if contig not encountered before
                    sample_singleton_pos[sample][contig] = []    #add a list for contigs for each sample
                sample_singleton_pos[sample][contig].append(pos)    #positions of singleton SNPs added to sample dict

            
total = 0     #counting all singletons stored, starting at 0 (quick check)
for sample in sample_singleton_pos.values():    #counting the values stored in nested dict
    for contig, positions in sample.items():    #grabbing contig name + position of each singleton mutation
        total += len(positions)               #counting how many positions are in list i.e. how many singletons present overall 
print("Total singleton positions stored:", total)  #e.g. for ST = 3250

#filtering for recombination
recombinant_pos = set()   #empty set storing recombinant positions
for sample in sample_singleton_pos:     #for each singleton stored 
    for contig, positions in sample_singleton_pos[sample].items():  #looping through contigs & list of singleton positions
        positions.sort()       #sorting the positions from smallest to largest, getting ready to run 2000bp window
        for i in range(len(positions)):   #looping through each position & using as a starting point 
            window = []   #creates list which will store recombinant positions if three or more present within 2000bp window
            for j in range(i + 1, len(positions)):  #starting from pos after i, looking at every next pos
                if positions[j] - positions[i] <= 2000:  #calculating distance between pos i and j
                    window.append(positions[j])   #if distance within 2000bp, add to window list
                else:
                    break   #if pos more than 2000bp, stop checking for SNPs 
            if len(window) >= 3:    #if 3/more pos (aka SNPs) within 2000bp, flag as recomb
                recombinant_pos.add((contig, positions[i]))   #add starting pos to recomb set
                for pos in window:      #loop through nearby pos within 2000bp window
                    recombinant_pos.add((contig, pos))   #add pos to recomb set

print("Total recombinant positions:", len(recombinant_pos))   #checking number of recomb positions/SNPs found 
print("Number of contigs with recombinant positions:", len(set(contig for contig, pos in recombinant_pos)))  #checking how many contigs have recomb present

with open("ST_recombinant_snps.csv", "w") as f:  #creating csv to store recombinant SNPs
    f.write("CHROM,POS\n")  #storing contig & position at which recomb flagged
    for contig, pos in recombinant_pos:   #going through all recomb-flagged pos
        f.write(f"{contig}, {pos}\n")    #adding contig & pos of recomb SNPs to csv
print("Recombinant SNPs saved to ST_recombinant_snps.csv") 

with open("ST_clean_singletons_test.csv", "w") as f:  #creating csv to store recomb-free singletons
    f.write("CHROM,POS\n")     #storing contig & pos 
    for sample in sample_singleton_pos:     #looping through genomes
        for contig, positions in sample_singleton_pos[sample].items():   #looping through contigs + positions within
            for pos in positions:    #looping through each singleton pos 
                if (contig, pos) not in recombinant_pos:   #if not flagged as recomb
                  f.write(f"{contig}, {pos}\n")   #add to clean singletons csv
print("Clean singletons saved to ST_clean_singletons.csv")


#creating 3 new vcfs to store all singletons, recomb singletons and clean singletons

with open("ST_raw.vcf", "r") as infile,   #using raw vcf as input 
     open("ST_singletons_only.vcf", "w") as singletons_out,  #generating vcf with all singletons
     open("ST_recombinant_snps.vcf", "w") as recomb_out,   #generating vcf with all recombinant singletons
     open("ST_clean_singletons.vcf", "w") as clean_out:      #generating vcf with non-recombinant singletons - will be using for analysis 

    for line in infile:  #going through raw vcf file
        if line.startswith("#"):   #copying headers into all three new vcf files 
            singletons_out.write(line)
            recomb_out.write(line)
            clean_out.write(line)
            continue
        entries = line.strip().split("\t")  #splitting lines into columns
        if len(entries) < 10:   #skip if less than 10 columns 
            continue
        if entries[6] != "PASS":   #filtering for "PASS" 
            continue
        ref = entries[3]  
        alt = entries[4]
        if len(ref) != 1 or len(alt) != 1:  #skipping multi-allelic sites
            continue
        if ref not in only_bases or alt not in only_bases:  #skipping rows if not include either of the valid bases
            continue

        sample_genotypes = entries[10:]   #picking out sample genotype cols 
        carriers = sum(1 for g in sample_genotypes if g == "1")  #counting no of samples carrying mut
        if carriers != 1:   #skipping if not singleton
            continue
        contig = entries[0]    #contig names
        pos = int(entries[1])  #position

        singletons_out.write(line)   #writing every singleton picked up into singletons_only vcf file

        if (contig,pos) in recombinant_pos:   #if position flagged as recomb
            recomb_out.write(line)   #write out recomb singletons in recombinant_snps vcf
        else:
            clean_out.write(line)   #everything else added to clean singletons vcf

print("VCF files saved!")


