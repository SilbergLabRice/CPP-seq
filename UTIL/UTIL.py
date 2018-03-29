import re, csv, numpy, math, scipy
import scipy.optimize
import scipy.stats as ss
import Bio.SeqIO

#Utilities for analyzing Transposon Insertions in NGS data set 
def TransposonGeneJunctionExtraction(filename, five_fw_tag, five_rev_tag, three_fw_tag, three_rev_tag, kmer = 54):
    # from Bio import SeqIO
    print('Processing dataset: ' + filename)
    print('Finding Transposon Tagged sequences...')
    
    #Extractions
    five_fw_tag = re.findall('([A,C,G,T]{' + str(kmer) + '})', five_fw_tag[::-1])[0][::-1]        	#grab last K-merbp adjacent to start codon 
    five_rev_tag = re.findall('([A,C,G,T]{' + str(kmer) + '})', five_rev_tag)[0]
    three_fw_tag = re.findall('([A,C,G,T]{' + str(kmer) + '})', three_fw_tag)[0]
    three_rev_tag = re.findall('([A,C,G,T]{' + str(kmer) + '})', three_rev_tag[::-1])[0][::-1] 		#grab last K-merbp adjacent to stop codon 

    print('\t' + 'Forward tag: ' + five_fw_tag)
    print('\t' + 'Reverse Complement tag: ' + five_rev_tag)
    print('\t' + 'Forward tag: ' + three_fw_tag)
    print('\t' + 'Reverse Complement tag: ' + three_rev_tag)

    #Number of records with 5' tag sequence (Start codon and 18aa tag)
    tagAFwdlist = []
    tagARevlist = []
    for record in Bio.SeqIO.parse(filename, "fastq"): 
        if re.search(five_fw_tag, str(record.seq)):
            tagAFwdlist.append(record.id)
        if re.search(five_rev_tag , str(record.seq)):
            tagARevlist.append(record.id)
    totalA = len(tagAFwdlist) + len(tagARevlist)
    print('\t' + "There are " + str(totalA) + ' records containing the 5 prime tag')

    #Number of records with 3' tag sequence (the stop codon and NotI tag)

    tagBFwdlist = []
    tagBRevlist = []
    for record in Bio.SeqIO.parse(filename, "fastq"):
        if re.search(three_fw_tag, str(record.seq)):
            tagBFwdlist.append(record.id)
        if re.search(three_rev_tag, str(record.seq)):
            tagBRevlist.append(record.id)
    totalB = len(tagBFwdlist) + len(tagBRevlist)        
    print('\t' + "There are " + str(totalB) + ' records containing the 3 prime tag')

    #Generate FASTQ with only tagged reads 
    wanted_ids = tagAFwdlist + tagARevlist + tagBFwdlist + tagBRevlist
    wanted_ids = list(set(wanted_ids)) #cleanup any duplicates 
    input_filename = filename 
    output_filename = 'Outputs/' + re.findall('\w.*/(\w.*).fastq', filename)[0] + '_tagged.fastq'
    
    fasta_index = Bio.SeqIO.index(input_filename, "fastq")
    count = 0
    total = len(fasta_index)
    output_handle = open(output_filename, "w")
    for identifier in wanted_ids:
        record = fasta_index[identifier]
        Bio.SeqIO.write(record, output_handle, "fastq")
        count = count + 1
    output_handle.close()
    
    print('\t' + str(count) + " records selected out of " + str(total))
    print('\t' +  str((float(count)/total)*100) + "% of the reads are useful for identifying insertion sites")
    print('Tagged Sequences written to FASTQ')
    print('\t' + output_filename)
    print('\n')

    return output_filename

def AdjacentMerExtraction(filename, five_fw_tag, five_rev_tag, three_fw_tag, three_rev_tag, kmer = 54):
    filename = 'Outputs/' + re.findall('\w.*/(\w.*).fastq', filename)[0] + '_tagged.fastq'
    
    five_fw_tag = re.findall('([A,C,G,T]{' + str(kmer) + '})', five_fw_tag[::-1])[0][::-1]        	#grab last K-merbp adjacent to start codon 
    five_rev_tag = re.findall('([A,C,G,T]{' + str(kmer) + '})', five_rev_tag)[0]
    three_fw_tag = re.findall('([A,C,G,T]{' + str(kmer) + '})', three_fw_tag)[0]
    three_rev_tag = re.findall('([A,C,G,T]{' + str(kmer) + '})', three_rev_tag[::-1])[0][::-1] 		#grab last K-merbp adjacent to stop codon
    
    print('Extracting Adjacent 5- & 11-mers...')
    tagAFwdlist = []
    tagARevlist = []
    tagAFwdinsert = []
    tagARevinsert = []
    tagBFwdlist = []
    tagBRevlist = []
    tagBFwdinsert = []
    tagBRevinsert = []
    orient1 = []
    orient2= []
    orient3= []
    orient4= []
    ID = []
    elevenmer = [] #verified to be unique of a sequence this length
    fivemer = []
    fivemerA = []
    fivemerB = []
    fivemerC = []
    fivemerD = []
    x = []
    y=[]

    for record in Bio.SeqIO.parse(filename, "fastq"):
        if re.search(five_fw_tag + '([A,C,G,T]{11})', str(record.seq)): 
            x = re.findall(five_fw_tag + '([A,C,G,T]{11})', str(record.seq))
            y =  re.findall(five_fw_tag + '([A,C,G,T]{5})', str(record.seq))
            tagAFwdlist.append(record.id)
            tagAFwdinsert.append(x[0])
            fivemerA.append(y[0])
            orient1.append('AF')
        if re.search('([A,C,G,T]{11})' + five_rev_tag, str(record.seq)):
            x = re.findall('([A,C,G,T]{11})' + five_rev_tag, str(record.seq))
            y =  re.findall('([A,C,G,T]{5})' + five_rev_tag, str(record.seq))
            tagARevlist.append(record.id)
            tagARevinsert.append(x[0])
            fivemerB.append(y[0])
            orient2.append('AR')
        if re.search('([A,C,G,T]{11})' + three_fw_tag, str(record.seq)):
            x = re.findall('([A,C,G,T]{11})' + three_fw_tag, str(record.seq))
            y =  re.findall('([A,C,G,T]{5})' + three_fw_tag, str(record.seq))
            tagBFwdlist.append(record.id)
            tagBFwdinsert.append(x[0])
            fivemerC.append(y[0])
            orient3.append('BF')
        if re.search(three_rev_tag + '([A,C,G,T]{11})', str(record.seq)): 
            x = re.findall(three_rev_tag + '([A,C,G,T]{11})', str(record.seq))
            y =  re.findall(three_rev_tag + '([A,C,G,T]{5})', str(record.seq))
            tagBRevlist.append(record.id)
            tagBRevinsert.append(x[0])
            fivemerD.append(y[0])
            orient4.append('BR')

    
    ID = tagAFwdlist + tagARevlist + tagBFwdlist + tagBRevlist
    elevenmer = tagAFwdinsert + tagARevinsert + tagBFwdinsert + tagBRevinsert
    orientation = orient1 + orient2 + orient3 + orient4
    fivemer = fivemerA + fivemerB + fivemerC + fivemerD
    print('\t' + str(len(elevenmer)) + " entries with gene matching 11mers flanking tag printed to CSV")
    print('\t' + str(len(numpy.unique(fivemer)))+ ' unique 5mers identified')

    rows = zip(ID, elevenmer, fivemer, orientation)
    output_filename = 'Outputs/' + re.findall('\w.*/(\w.*)_tagged.fastq', filename)[0] + '_raw.csv'
    with open(output_filename, 'wb') as f:
        writer = csv.writer(f)
        writer.writerow(['ID', '11mer','5mer', 'Tag'])
        for row in rows:
            writer.writerow(row)

    print('Written to CSV')
    print('\t' + output_filename)
    print('\n')

    return output_filename

def PositionExtraction(filename, geneName, gene_fwd, gene_rev):
    print('Identifying Transposon Insertion Positions...')
    
    filename = 'Outputs/' + re.findall('\w.*/(\w.*).fastq', filename)[0] + '_raw.csv'
    
    possibleSites = numpy.genfromtxt('Outputs/' + geneName + '_PossibleSites.csv', delimiter = ',', names=True, dtype=None)
    possibleIndex = []
    possibleElevenmer = []
    possibleOligo = []

    for i in possibleSites:
        possibleIndex.append(i[0])
        possibleOligo.append(i[1])
        possibleElevenmer.append(i[3])


    data = numpy.genfromtxt(filename, delimiter = ',', names=True, dtype=None)
    tag1 = data['11mer'][data['Tag'] == 'AF'] #Fwd is the only way to make a real protein
    tag2 = data['11mer'][data['Tag'] == 'AR'] #Revcomp is the only way to make a real protein
    tag3 = data['11mer'][data['Tag'] == 'BF'] #Fwd is the only way to make a real protein
    tag4 = data['11mer'][data['Tag'] == 'BR'] #Revcomp is the only way to make a real protein


    frame = []
    residue = []
    bp = []
    direction = []
    orientation = []
    oligo = []


    for eleven in tag1:

        #search for index of matching elevenmer in FWD/REVCOMP look up table
        nucleotide = [possibleIndex[ind] for ind, x in enumerate(possibleElevenmer) if x == eleven]
        fragment = [possibleOligo[ind] for ind, x in enumerate(possibleElevenmer) if x == eleven]

        #determine if there is 0,1, or multiple non-unique matches
        if len(nucleotide) == 0:
            #print('No Match Found')
            frame.append(0)
            residue.append(0)
            bp.append(0)
            oligo.append(0)
            direction.append('No Match')
            orientation.append('No Match')
            continue
        if len(nucleotide) > 1:
            print('Multiple Matches Found')
            frame.append(0)
            residue.append(0)
            bp.append(0)
            oligo.append(0)
            direction.append('Multiple Matches')
            orientation.append('Multiple Matches')
            continue 
        nuc = nucleotide[0]  #used for orientation 
        nucleotide = nucleotide[0] #used for index adjustment for other tags
        oligo.append(fragment[0])

        #search Fwd look-up table
        if nuc <= len(gene_fwd):
            nucleo=nucleotide
            if nucleo < 0:
                nucleo =  -nucleo
            if nucleo > len(gene_fwd):
                nucleo =   nucleo - len(gene_fwd) 
            fm = (nucleo) % 3
            if fm == 0:
                fm = 3

            res = int(math.ceil(float(nucleo)/3))

            frame.append(fm)    
            residue.append(res)
            bp.append(nucleo)
            direction.append('F')
            orientation.append('Fwd')
            continue
        
        #search Revcomp look-up table
        if nuc > len(gene_fwd):
            nucleo = nucleotide - len(gene_fwd) + 10 -4 #adjust index for flipped insertion 
            if nucleo < 0:
                nucleo = -nucleo 
            if nucleo > len(gene_fwd):
                nucleo = nucleo - len(gene_fwd) 

            fm = (nucleo) % 3
            if fm == 0:
                fm = 3

            res = int(math.ceil(float(nucleo)/3))

            frame.append(fm)
            residue.append(res)
            bp.append(nucleo)
            direction.append('NF')
            orientation.append('Rev')
            continue
    
    for eleven in tag2:

        #search for index of matching elevenmer in FWD/REVCOMP look up table
        nucleotide = [possibleIndex[ind] for ind, x in enumerate(possibleElevenmer) if x == eleven]
        fragment = [possibleOligo[ind] for ind, x in enumerate(possibleElevenmer) if x == eleven]

        #determine if there is 0,1, or multiple non-unique matches
        if len(nucleotide) == 0:
            #print('No Match Found')
            frame.append(0)
            residue.append(0)
            bp.append(0)
            oligo.append(0)
            direction.append('No Match')
            orientation.append('No Match')
            continue
        if len(nucleotide) > 1:
            print('Multiple Matches Found')
            frame.append(0)
            residue.append(0)
            bp.append(0)
            oligo.append(0)
            direction.append('Multiple Matches')
            orientation.append('Multiple Matches')
            continue 
        
        nuc = nucleotide[0]
        nucleotide = nucleotide[0]-10
        oligo.append(fragment[0])



        #search Fwd look-up table
        if nuc <= len(gene_fwd):
            nucleo = nucleotide + 10 + 6

            if nucleo < 0:
                nucleo =  -nucleo
            if nucleo > len(gene_fwd):
                nucleo =   nucleo - len(gene_fwd)

            fm = (nucleo) % 3
            if fm == 0:
                fm = 3

            res = int(math.ceil(float(nucleo)/3))

            frame.append(fm)    
            residue.append(res)
            bp.append(nucleo)
            direction.append('NF')
            orientation.append('Fwd')
            continue
            
        #search Revcomp look-up table
        if nuc > len(gene_fwd):
            nucleo = nucleotide - len(gene_fwd) + 10 #adjust index for flipped insertion 

            if nucleo < 0:
                nucleo = -nucleo 
            if nucleo > len(gene_fwd):
                nucleo =  nucleo - len(gene_fwd)

            fm = (nucleo) % 3
            if fm == 0:
                fm = 3

            res = int(math.ceil(float(nucleo)/3))

            frame.append(fm)
            residue.append(res)
            bp.append(nucleo)
            direction.append('F')
            orientation.append('Rev')
            continue
            
    for eleven in tag3:

        #search for index of matching elevenmer in FWD/REVCOMP look up table
        nucleotide = [possibleIndex[ind] for ind, x in enumerate(possibleElevenmer) if x == eleven]
        fragment = [possibleOligo[ind] for ind, x in enumerate(possibleElevenmer) if x == eleven]

        #determine if there is 0,1, or multiple non-unique matches
        if len(nucleotide) == 0:
            #print('No Match Found')
            frame.append(0)
            residue.append(0)
            bp.append(0)
            oligo.append(0)
            direction.append('No Match')
            orientation.append('No Match')
            continue
        if len(nucleotide) > 1:
            print('Multiple Matches Found')
            frame.append(0)
            residue.append(0)
            bp.append(0)
            oligo.append(0)
            direction.append('Multiple Matches')
            orientation.append('Multiple Matches')
            continue 

        nuc = nucleotide[0]
        nucleotide = nucleotide[0] + 6
        oligo.append(fragment[0])


        #search Fwd look-up table
        if nuc <= len(gene_fwd):
            nucleo = nucleotide
            if nucleo < 0:
                nucleo =  6 + nucleo
            if nucleo > len(gene_fwd):
                nucleo =   6 + nucleo - len(gene_fwd)
            fm = (nucleo) % 3
            if fm == 0:
                fm = 3

            res = int(math.ceil(float(nucleo)/3))

            frame.append(fm)    
            residue.append(res)
            bp.append(nucleo)
            direction.append('F')
            orientation.append('Fwd')
            continue
            
        #search Revcomp look-up table
        if nuc > len(gene_fwd):
            nucleo = nucleotide - len(gene_fwd) + 10 - 10 - 6 #adjust index for flipped insertion 

            if nucleo < 0:
                nucleo = 6 + nucleo 
            if nucleo > len(gene_fwd):
                nucleo = 6 + nucleo - len(gene_fwd) 

            fm = (nucleo) % 3
            if fm == 0:
                fm = 3

            res = int(math.ceil(float(nucleo)/3))

            frame.append(fm)
            residue.append(res)
            bp.append(nucleo)
            direction.append('NF')
            orientation.append('Rev')
            continue
            
    for eleven in tag4:

        #search for index of matching elevenmer in FWD/REVCOMP look up table
        nucleotide = [possibleIndex[ind] for ind, x in enumerate(possibleElevenmer) if x == eleven]
        fragment = [possibleOligo[ind] for ind, x in enumerate(possibleElevenmer) if x == eleven]

        #determine if there is 0,1, or multiple non-unique matches
        if len(nucleotide) == 0:
            #print('No Match Found')
            frame.append(0)
            residue.append(0)
            bp.append(0)
            oligo.append(0)
            direction.append('No Match')
            orientation.append('No Match')
            continue
        if len(nucleotide) > 1:
            print('Multiple Matches Found')
            frame.append(0)
            residue.append(0)
            bp.append(0)
            oligo.append(0)
            direction.append('Multiple Matches')
            orientation.append('Multiple Matches')
            continue 

        nuc = nucleotide[0]
        nucleotide = nucleotide[0] - 4
        oligo.append(fragment[0])

        #search Fwd look-up table
        if nuc <= len(gene_fwd):
            nucleo = nucleotide + 4

            if nucleo < 0:
                nucleo =  6 + nucleo
            if nucleo > len(gene_fwd):
                nucleo =   6 + nucleo - len(gene_fwd)

            fm = (nucleo) % 3
            if fm == 0:
                fm = 3

            res = int(math.ceil(float(nucleo)/3))

            frame.append(fm)    
            residue.append(res)
            bp.append(nucleo)
            direction.append('NF')
            orientation.append('Fwd')
            continue
            
        #search Revcomp look-up table
        if nuc > len(gene_fwd):
            nucleo = nucleotide - len(gene_fwd) + 10 #adjust index for flipped insertion 

            if nucleo < 0:
                nucleo =  6 + nucleo
            if nucleo > len(gene_fwd):
                nucleo =   6 + nucleo - len(gene_fwd) 

            fm = (nucleo) % 3
            if fm == 0:
                fm = 3

            res = int(math.ceil(float(nucleo)/3))

            frame.append(fm)
            residue.append(res)
            bp.append(nucleo)
            direction.append('F')
            orientation.append('Rev')
            continue
            
    rows = zip(data['ID'], oligo, data['11mer'], data['5mer'], data['Tag'], bp, frame, residue, direction, orientation)
    output_filename = 'Outputs/' + re.findall('\w.*/(\w.*)_raw.csv', filename)[0] + '_positions.csv'
    with open(output_filename, 'wb') as f:
        writer = csv.writer(f)
        writer.writerow(['ID', 'oligo', '11mer', '5mer', 'Tag', 'position', 'frame', 'residue', 'direction', 'gene_strand'])
        for row in rows:
            writer.writerow(row)

    res = numpy.array(residue)
    sample = res[res != 0]

    b = numpy.array(direction)
    F = b<'NF'
    NF = b>'F'

    forward = res[F]
    reverse = res[NF]

    print('\t' + str(len(res)) + ' total reads') #total residues
    print('\t' + str(len(forward)) + ' forward reads')
    print('\t' + str(len(reverse))  + ' reverse reads' ) #subselection

    print('Written to CSV')
    print('\t' + output_filename)
    print('\n')
    
    return output_filename

def NucleotideInsertionSiteCounter(file, gene_fwd): 
    print('Counting Transposon Insertion Positions...')
    filename = 'Outputs/' + re.findall('\w.*/(\w.*).fastq', file)[0] + '_positions.csv'
    
    possibleSites = numpy.genfromtxt('Outputs/' + re.findall('\w.*/(\w.*)_\w.*.fastq', file)[0] + '_PossibleSites.csv', delimiter = ',', names=True, dtype=None)
    possibleOligo = []
    possibleFivemer = []
    possibleElevenmer = []
    

    for i in possibleSites:
        possibleOligo.append(i[1])
        possibleFivemer.append(i[2])
        possibleElevenmer.append(i[3])
    
    data = numpy.genfromtxt(filename, delimiter = ',', names=True, dtype=None)
    functional = data[data['direction'] == 'F']
    nonfunctional = data[data['direction'] == 'NF']
    
    print('\t There are ' + str(len(functional)) + ' parallel insertions')
    print('\t There are ' + str(len(nonfunctional)) + ' antiparallel insertions')

    index = []
    oligo = []
    eleven = []
    five = []
    counts = []
    
    for i in range(0, len(gene_fwd)):
        index.append(i + 1)
        oligo.append(possibleOligo[i])
        eleven.append(possibleFivemer[i])
        five.append(possibleElevenmer[i])
        counts.append(sum(functional['position'] == i + 1))
        
    for i in range(0, len(gene_fwd)):
        i2 = i + len(gene_fwd)
        index.append(i2 + 1)
        oligo.append(possibleOligo[i2])
        eleven.append(possibleFivemer[i2])
        five.append(possibleElevenmer[i2])
        counts.append(sum(nonfunctional['position'] == i + 1))


    output_filename = 'Outputs/' + re.findall('\w.*/(\w.*)_positions.csv', filename)[0] + '_NucleotidePositionCounts.csv'
    rows = zip(index, oligo, eleven, five, counts)
    with open(output_filename, 'wb') as f:
        writer = csv.writer(f)
        writer.writerow(['index', 'oligo', '5mer', '11mer', 'Parallel_Count', 'Antiparallel_Count'])
        for row in rows:
            writer.writerow(row)
            
    print('Written to CSV')
    print('\t' + output_filename)
    print('\n')
    
def AminoAcidStartSiteCounter(file, gene_fwd):
    print('Counting new N-termini...')
    filename = 'Outputs/' + re.findall('\w.*/(\w.*).fastq', file)[0] + '_positions.csv'
    
    data = numpy.genfromtxt(filename, delimiter = ',', names=True, dtype=None)
    functional = data[data['direction'] == 'F']
    nonfunctional = data[data['direction'] == 'NF']
    nomatch = data[data['direction'] == 'No Match'] 
    multiple = data[data['direction'] == 'Multiple Matches']

    functional_inframe = functional[functional['frame'] == 1]
    nonfunctional_inframe = nonfunctional[nonfunctional['frame'] == 1]

    functional_frametwo = functional[functional['frame'] == 2]
    nonfunctional_frametwo = nonfunctional[nonfunctional['frame'] == 2]

    functional_framethree = functional[functional['frame'] == 3]
    nonfunctional_framethree = nonfunctional[nonfunctional['frame'] == 3]

    print('\t There are ' + str(len(functional_inframe)) + ' parallel in frame insertions')
    print('\t There are ' + str(len(nonfunctional_inframe)) + ' antiparallel in frame insertions')
    print('\t There are ' + str(len(functional_frametwo)) + ' parallel +1 frame insertions')
    print('\t There are ' + str(len(nonfunctional_frametwo)) + ' antiparallel +1 frame insertions')
    print('\t There are ' + str(len(functional_framethree)) + ' parallel -1 frame insertions')
    print('\t There are ' + str(len(nonfunctional_framethree)) + ' antiparallel -1 frame insertions')
    print('\t There are ' + str(len(nomatch)) + ' insertions that do not match VARIABLE region')
    print('\t There are ' + str(len(multiple)) + ' insertions that with multiple matches VARIABLE region')

    aa = numpy.zeros((len(gene_fwd)/3,7))
    for i in range(0, len(aa)): #len(aa)
        aa[i,0] = i + 1
        aa[i,1] = sum(functional_inframe['residue'] == i + 1)
        aa[i,2] = sum(functional_frametwo['residue'] == i + 1)
        aa[i,3] = sum(functional_framethree['residue'] == i + 1)
        aa[i,4] = sum(nonfunctional_inframe['residue'] == i + 1)
        aa[i,5] = sum(nonfunctional_frametwo['residue'] == i + 1)
        aa[i,6] = sum(nonfunctional_framethree['residue'] == i + 1)

    print(len(aa))


    output_filename = 'Outputs/' + re.findall('\w.*/(\w.*)_positions.csv', filename)[0] + '_ProteinNTerminiCounts.csv'
    rows = zip(aa[:,0], aa[:,1], aa[:,2], aa[:,3], aa[:,4], aa[:,5], aa[:,6])
    with open(output_filename, 'wb') as f:
        writer = csv.writer(f)
        writer.writerow(['index', 'Parallel_F1_Count', 'Parallel_F2_Count', 'Parallel_F3_Count', 'Antiparallel_F1_Count', 'Antiparallel_F2_Count', 'Antiparallel_F3_Count'])
        for row in rows:
            writer.writerow(row)
            
    print('Written to CSV')
    print('\t' + output_filename)
    print('\n')


def FoldChangeCalculator(geneName):
    print('Calculating fold changes for each frame of ' + geneName + '...')
    #Takes AAStartSites for each frame and calculates log2(fold) change for selected/unselected 
    selected_file = 'Outputs/' + geneName + '_selected_ProteinNTerminiCounts.csv'
    unselected_file = 'Outputs/' + geneName + '_unselected_ProteinNTerminiCounts.csv'

    globals()[re.findall('Outputs/\w.*_(\w.*)_\w.*.csv', selected_file)[0]] = numpy.genfromtxt(selected_file, delimiter = ',', names=True, dtype=None)
    globals()[re.findall('Outputs/\w.*_(\w.*)_\w.*.csv', unselected_file)[0]] = numpy.genfromtxt(unselected_file, delimiter = ',', names=True, dtype=None)

    index = selected['index']

    for i in ['1', '2', '3']:
        selected_Parallel_frame = selected['Parallel_F' + i + '_Count']
        selected_Antiparallel_frame = selected['Antiparallel_F' + i + '_Count']

        unselected_Parallel_frame = unselected['Parallel_F' + i + '_Count']
        unselected_Antiparallel_frame = unselected['Antiparallel_F' + i + '_Count']

        unobserved_selected_Parallel = [ind for ind, x in enumerate(selected_Parallel_frame) if x == 0]
        unobserved_unselected_Parallel = [ind for ind, x in enumerate(unselected_Parallel_frame) if x == 0]

        unobserved_selected_Antiparallel = [ind for ind, x in enumerate(selected_Antiparallel_frame) if x == 0]
        unobserved_unselected_Antiparallel = [ind for ind, x in enumerate(unselected_Antiparallel_frame) if x == 0]

        noLongerObserved_Parallel = numpy.setdiff1d(unobserved_selected_Parallel, unobserved_unselected_Parallel)
        noLongerObserved_Antiparallel = numpy.setdiff1d(unobserved_selected_Antiparallel, unobserved_unselected_Antiparallel)

        #Eliminate Zeros that throw errors
        unselected_Parallel_frame = [1 if x==0 else x for x in unselected_Parallel_frame]
        unselected_Antiparallel_frame = [1 if x==0 else x for x in unselected_Antiparallel_frame]

        print('\n')
        print('\t' +'Are there any residues observed in the Selected but not the Unselected library? - Parallel: ')
        print(len(unobserved_unselected_Parallel ) != len(numpy.intersect1d(unobserved_selected_Parallel, unobserved_unselected_Parallel)))
        print('\t' +'Are there any residues observed in the Selected but not the Unselected library? - Antiparallel: ')
        print(len(unobserved_unselected_Antiparallel ) != len(numpy.intersect1d(unobserved_selected_Antiparallel, unobserved_unselected_Antiparallel)))

        print('\n')
        print('\t' +'There are ' + str(len(noLongerObserved_Parallel)) + ' residues no longer observed in the selected library - Parallel')
        print('\t' +'There are ' + str(len(noLongerObserved_Antiparallel)) + ' residues no longer observed in the selected library - Antiparallel')

        globals()['F' + i + '_Parallel_fold'] = selected_Parallel_frame/unselected_Parallel_frame
        globals()['F' + i + '_Antiparallel_fold'] = selected_Antiparallel_frame/unselected_Antiparallel_frame

        #Eliminate Zeros that throw errors
        globals()['F' + i + '_Parallel_fold'] = [0.5 if x==0 else x for x in globals()['F' + i + '_Parallel_fold']]
        globals()['F' + i + '_Antiparallel_fold'] = [0.5 if x==0 else x for x in globals()['F' + i + '_Antiparallel_fold']]

        #calculate log2(fold change)
        globals()['F' + i + '_Parallel_fold_log'] = numpy.log2(globals()['F' + i + '_Parallel_fold'])
        globals()['F' + i + '_Antiparallel_fold_log'] = numpy.log2(globals()['F' + i + '_Antiparallel_fold'])

        #replace former zeros with 'nan' to mark no longer observed positions that are infinitely diluted
        globals()['F' + i + '_Parallel_fold_log'][noLongerObserved_Parallel] = 'nan' 
        globals()['F' + i + '_Antiparallel_fold_log'][noLongerObserved_Antiparallel] = 'nan'

    output_filename = 'Outputs/' + geneName + '_FoldChange.csv'
    rows = zip(index, F1_Parallel_fold_log, F2_Parallel_fold_log, F3_Parallel_fold_log, F1_Antiparallel_fold_log, F2_Antiparallel_fold_log, F3_Antiparallel_fold_log)
    with open(output_filename, 'wb') as f:
        writer = csv.writer(f)
        writer.writerow(['index', 'Parallel_F1_Fold', 'Parallel_F2_Fold', 'Parallel_F3_Fold', 'Antiparallel_F1_Fold', 'Antiparallel_F2_Fold', 'Antiparallel_F3_Fold'])
        for row in rows:
            writer.writerow(row)

    print('Written to CSV')
    print('\t' + output_filename)
    print('\n')
    
    return rows 

def FisherTest(geneName):
    geneName = geneName
    print('Calculating Fisher Exact Test on ' + geneName + '...')

    selected_file = 'Outputs/' + geneName + '_selected_ProteinNTerminiCounts.csv'
    unselected_file = 'Outputs/' + geneName + '_unselected_ProteinNTerminiCounts.csv'

    globals()[re.findall('Outputs/\w.*_(\w.*)_\w.*.csv', selected_file)[0]] = numpy.genfromtxt(selected_file, delimiter = ',', names=True, dtype=None)
    globals()[re.findall('Outputs/\w.*_(\w.*)_\w.*.csv', unselected_file)[0]] = numpy.genfromtxt(unselected_file, delimiter = ',', names=True, dtype=None)

    for k in [1,2,3]:
        pv = []

        for i in range(0,len(selected)):
            selected_P = selected[i][k]
            selected_AP = selected[i][k+3]
            unselected_P = unselected[i][k]
            unselected_AP = unselected[i][k+3]            
            
            selected_sample = [selected_P, selected_AP]
            unselected_sample = [unselected_P, unselected_AP]
            sample = [unselected_sample, selected_sample]
            
            sample = [1 if x == 0 else x for x in sample]
            
            ft = ss.fisher_exact(sample)[1]
            if ft == 0:
                ft = 10**-300
 
            pv.append(numpy.log10(ft))

        globals()['F' + str(k) + '_pv'] = pv
 
    index = selected['index']
    output_filename = 'Outputs/' + geneName + '_FisherTest.csv'
    #note fishers test the P-values are the same for both P/AP since contingency table
    rows = zip(index, F1_pv, F2_pv, F3_pv, F1_pv, F2_pv, F3_pv)
    with open(output_filename, 'wb') as f:
        writer = csv.writer(f)
        writer.writerow(['index', 'Frame1 P - log10(pvalue)', 'Frame2 P - log10(pvalue)', 'Frame3 P - log10(pvalue)', 'Frame1 AP - log10(pvalue)', 'Frame2 AP - log10(pvalue)', 'Frame3 AP - log10(pvalue)'])
        for row in rows:
            writer.writerow(row)

    print('Written to CSV')
    print('\t' + output_filename)
    print('\n')
    
def NBTest(geneName):    
    geneName = geneName
    print('Calculating Negative Binomial Test on ' + geneName + '...')

    selected_file = 'Outputs/' + geneName + '_selected_ProteinNTerminiCounts.csv'
    unselected_file = 'Outputs/' + geneName + '_unselected_ProteinNTerminiCounts.csv'


    #import data
    globals()[re.findall('\w.*_(\w.*)_\w.*.csv', selected_file)[0]] = numpy.genfromtxt(selected_file, delimiter = ',', names=True, dtype=None)
    globals()[re.findall('\w.*_(\w.*)_\w.*.csv', unselected_file)[0]] = numpy.genfromtxt(unselected_file, delimiter = ',', names=True, dtype=None)

    #loop over reading frames
    for k in [1, 2, 3]:

        #split data
        unselected_index = []
        US_P = []
        US_AP = []
        S_P = []
        S_AP = []


        for i in range(0, len(selected)):
            unselected_index.append(unselected[i][0])
            S_P.append(selected[i][k])
            S_AP.append(selected[i][k+3])
            US_P.append(unselected[i][k])
            US_AP.append(unselected[i][k+3])

        ################################
        #####Parameter Estimation#######
        ################################

        #classify and remove low count positions (<20 reads) from unselected-AP for estimation 
        remove_unobserved = True
        observation_threshold = 20

        if remove_unobserved == False:
            US_AP_observed = range(0, len(unselected_index))

        if remove_unobserved == True:
            
            US_AP_observed = [int(unselected_index[index])-1 for index, value in enumerate(US_AP) if value > observation_threshold]
            US_AP_unobserved = [int(unselected_index[index])-1 for index, value in enumerate(US_AP) if value <= observation_threshold]

        print('Postions with >20 reads: ' + str(len(US_AP_observed)) + '/' + str(len(unselected_index)))


        ### Estimate Beta ###
        AP_ratio = [S_AP[i]/US_AP[i] for i in US_AP_observed] 

        beta = numpy.median(AP_ratio)
        print('Beta: ' + str(beta))

        ### Estimate Gamma ###
        def fit_gamma(gamma): 
            like = 0
            for i in range(0,len(US_AP_observed)):
                x = US_AP[US_AP_observed[i]]
                y = S_AP[US_AP_observed[i]] 

                def ll(m):
                    r1 = m/(gamma-1)
                    r2 = beta * r1
                    p = 1/gamma

                    l1 = scipy.stats.nbinom.logpmf(x, r1, p)
                    l2 = scipy.stats.nbinom.logpmf(y, r2, p)

                    return -(l1+l2)

                #initial guess, matching R optimize function 
                bnds = (1e-8, 1e6)
                opt = scipy.optimize.minimize_scalar(ll, bounds = bnds, method = 'bounded') 
                like = like + opt.fun

            return like

        #initial guess, matching R optimize function
        bnds = (1.01, 100)
        gamma = scipy.optimize.minimize_scalar(fit_gamma, bounds = bnds, method = 'bounded').x 

        ### Output the estimated values for beta and gamma
        print("Gamma = " + str(gamma))

        ############## The functions for estimating mean and variance of log(FC) ###########
        def mu_log(r, p, c):
            mu = p * r / (1-p)
            v = p * r / (1-p) / (1-p)
            ml = numpy.log(mu + c) - v / 2 / (mu + c) / (mu + c)
            return ml

        def var_log(r, p, c):
            mu = p * r / (1-p)
            v = p*  r / (1-p) / (1-p)
            m4 = p * r * (p * p+4 * p+3 * p * r+1) / ((1-p)**4)
            v_log  = v/(mu + c)/(mu + c) + (m4 - v * v)/4/( (mu + c)**4)
            return v_log

        def mu_fc(us, p, c):
            r1 = us/(gamma - 1)
            r2 = r1 * beta
            mu1 = mu_log(r1, p, c)
            mu2 = mu_log(r2, p, c)
            return mu2-mu1

        def var_fc(us, p, c):
            r1 = us/(gamma - 1)
            r2 = r1 * beta
            v1 = var_log(r1, p, c)
            v2 = var_log(r2, p, c)
            return v1+v2

        def est_m(x, y):

            if ([x, y] == [0, 0]):
                return 0

            def ll(m):
                r1 = m/(gamma-1)
                r2 = beta * r1
                p = 1/gamma
                l1 = scipy.stats.nbinom.logpmf(x, r1, p)
                l2 = scipy.stats.nbinom.logpmf(y, r2, p)

                return -(l1+l2)    
            
            bnds = (1e-8, 1e10)
            opt = scipy.optimize.minimize_scalar(ll, bounds = bnds, method = 'bounded')

            return opt.x

        ## Function for computing p-value using negative binomial distribution ##
        def compute_pv(x, y, m, beta, gamma):
            r1 = m/(gamma - 1)
            r2 = beta * r1
            p = 1/gamma
            l1 = scipy.stats.nbinom.logpmf(x, r1, p)
            l2 = scipy.stats.nbinom.logpmf(y, r2, p)
            ll = l1 + l2
            sum_ll = ll 
            sum_small = ll
            s = x + y 
            for j in range(0, int(s)):
                if (j == y):
                    continue
                l1 = scipy.stats.nbinom.logpmf(s-j, r1, p)
                l2 = scipy.stats.nbinom.logpmf(j, r2, p)
                l = l1 + l2
                sum_ll = numpy.logaddexp(sum_ll, l)
                if l < ll:
                    sum_small = numpy.logaddexp(sum_small, l)
            log10_pv = (sum_small - sum_ll)/numpy.log(10)

            return log10_pv 

        def BH_correction_log(pv, alpha):
            #Benjamani Hochberg correction for log10 pvalues
            #sort pvalues and store originial index
            pv_sort = numpy.sort(pv)
            pv_index = numpy.argsort(pv)

            #Calculate BH critical values
            BH_critical = [(float(i)/len(pv_sort))*alpha for i in range(1, len(pv_sort)+1)]

            #Count significant values 
            largest_index = []
            for i in range(0, len(pv_sort)):
                if pv_sort[i] <= numpy.log10(BH_critical[i]):
                    largest_index.append(i)

            if len(largest_index) >= 1:
                nsig = max(largest_index)+1
            else:
                nsig = 0

            print(str(nsig) + ' values are significant after BH correction with a FDR of: ' + str(alpha))    


            #Adjust pvalues and reindex to original index
            pv_adj = []
            for i in range(1, len(pv_sort)+1):
                if i < len(pv_sort):
                    if (pv_sort[i-1]+numpy.log10(float(len(pv_sort))/i)) <= (pv_sort[i]+numpy.log10(float(len(pv_sort))/(i+1))):
                        pv_adj.append(pv_sort[i-1]+numpy.log10(float(len(pv_sort))/i))
                    if (pv_sort[i-1]+numpy.log10(float(len(pv_sort))/i)) > (pv_sort[i]+numpy.log10(float(len(pv_sort))/(i+1))):
                        pv_adj.append((pv_sort[i]+numpy.log10(float(len(pv_sort))/(i+1))))
                if i == len(pv_sort):
                    pv_adj.append(pv_sort[i-1]+numpy.log10(float(len(pv_sort))/i))

            pv_adj_reindex = numpy.zeros(len(pv_index))
            for i in range(0, len(pv_index)):
                if float(pv[pv_index[i]]) == 0.0:
                    pv_adj_reindex[pv_index[i]] = 0.0
                if float(pv[pv_index[i]]) == 1.0:
                    pv_adj_reindex[pv_index[i]] = 1.0
                if float(pv[pv_index[i]]) != 0.0 and float(pv[pv_index[i]]) != 1.0: 
                    pv_adj_reindex[pv_index[i]] = pv_adj[i]

            return pv_adj_reindex

        ################################
        #### Mean & StDev Estimation ###
        ################################

        c = 1  ### Feel free to play with this parameter. It does influence the significance...
        p = 1 - 1/gamma

        #classify and remove low count positions (0 reads) from unselected-AP for estimation 
        observation_threshold = 0
        AP_sel_na = [int(unselected_index[index])-1 for index, value in enumerate(S_AP) if value == observation_threshold]
        AP_unsel_na = [int(unselected_index[index])-1 for index, value in enumerate(US_AP) if value == observation_threshold]

        na = list(set(AP_sel_na) & set(AP_unsel_na))

        #if len(na) > 0 :
        #    AP_sel = numpy.delete(S_AP, na)
        #    AP_unsel = numpy.delete(US_AP, na)
        
        AP_sel = S_AP
        AP_unsel = US_AP
        
        #adjust by constant 
        AP_sel_fc = [i + c for i in AP_sel]
        AP_unsel_fc = [i + c for i in AP_unsel]

        #Calculate Fold Change
        fc = [numpy.log2((AP_sel_fc[i]+c)/(AP_unsel_fc[i]+c)) for i in range(0, len(AP_sel))]

        #estimate mean/stdev at each abundance level
        m = []
        mean_fc = []
        sd_fc = []

        for i in range(0,len(AP_sel)):
            mean_est = est_m(AP_unsel[i], AP_sel[i])
            m.append(mean_est)
            mean_fc.append(mu_fc(mean_est, p, c))
            sd_fc.append(numpy.sqrt(var_fc(mean_est, p, c)))

        #log2 transform data for plotting
        mean_fc = [i/numpy.log(2) for i in mean_fc]
        sd_fc = [i/numpy.log(2) for i in sd_fc]

        #sort for plotting mean/sd lines
        data = numpy.array(zip(m, fc, mean_fc, sd_fc), dtype = [('m', 'float64'), ('fc', 'float64'), ('mean', 'float64'), ('sd', 'float64')])
        data_ordered = numpy.sort(data, order = 'm')



        ################################
        ##### Calculating P-values #####
        ################################

        default_c = 0.5 ## this parameter is not that important 

        #Copy Count data
        AP_sel = S_AP
        AP_unsel = US_AP
        P_sel = S_P
        P_unsel = US_P

        #Calculate Fold Change log2(Selected/Unselected)
        fc_ap = [numpy.log2((AP_sel[i]+c)/(AP_unsel[i]+c)) for i in range(0, len(AP_sel))]
        fc_p = [numpy.log2((P_sel[i]+c)/(P_unsel[i]+c)) for i in range(0, len(P_unsel))]

        #Containers for p values
        pv_P = []
        pv_AP = []

        #calculate p-values 
        for i in range(0, len(AP_sel)):
            p1 = 0
            p2 = 0

            #Antiparallel p-value
            if (AP_sel[i] + AP_unsel[i]) == 0:
                p1 = 1.0 #'NA'
                fc_ap[i] = -9

            else:
                fc1 = numpy.log((AP_sel[i] + c)/(AP_unsel[i] + c))
                m1 = est_m(AP_unsel[i], AP_sel[i])
                mu1 = mu_fc(m1, p, c)
                v1 = var_fc(m1, p, c)
                p1 = compute_pv(AP_unsel[i], AP_sel[i], m1, beta, gamma)

            #Parallel p-value
            if (P_sel[i] + P_unsel[i] == 0):
                p2 = 1.0 #'NA'
                fc_p[i] = -9

            else: 
                fc2 = numpy.log((P_sel[i] + c)/(P_unsel[i] + c))
                m2 = est_m(P_unsel[i], P_sel[i])
                mu2 = mu_fc(m2, p, c)
                v2 = var_fc(m2, p, c)
                p2 = compute_pv(P_unsel[i], P_sel[i], m2, beta, gamma)
            
            #Adjust no-longer oberserved positions
            if AP_sel[i] == 0:
                fc_ap[i] = -9
            if P_sel[i] == 0:
                fc_p[i] = -9

            pv_AP.append(p1)
            pv_P.append(p2)

        #Multiple Comparison Correction using Benjamani-Hochberg 
        sig_thresh = 0.01
        stats_corrected_parallel = BH_correction_log(pv_P, alpha = sig_thresh)
        stats_corrected_antiparallel = BH_correction_log(pv_AP, alpha = sig_thresh)

        #Count significant pvalues
        counter = 0
        sig_thresh = -2
        for i in stats_corrected_parallel:
            if i <= sig_thresh:
                counter = counter + 1

        print('Look here significant P pv: ' + str(counter))
        globals()['F' + str(k) + '_P_pv'] = stats_corrected_parallel
        globals()['F' + str(k) + '_AP_pv'] = stats_corrected_antiparallel
        globals()['F' + str(k) + '_abund'] = data_ordered['m']
        globals()['F' + str(k) + '_mean'] = data_ordered['mean']
        globals()['F' + str(k) + '_sd'] = data_ordered['sd']


    pdat = zip(unselected_index, F1_P_pv, F2_P_pv, F3_P_pv, F1_AP_pv, F2_AP_pv, F3_AP_pv)
    pdat.insert(0,('index', 'Frame1 P - log10(pvalue)', 'Frame2 P - log10(pvalue)', 'Frame3 P - log10(pvalue)', 'Frame1 AP - log10(pvalue)', 'Frame2 AP - log10(pvalue)', 'Frame3 AP - log10(pvalue)'))
    
    fit = zip(F1_abund, F1_mean, F1_sd, F2_abund, F2_mean, F2_sd, F3_abund, F3_mean, F3_sd)
    fit.insert(0,('Frame1 - Abundance', 'Frame1 - mean', 'Frame1 - sd', 'Frame2 - Abundance', 'Frame2 - mean', 'Frame2 - sd', 'Frame3 - Abundance', 'Frame3 - mean', 'Frame3 - sd'))
    
    #Print P-vales to txt 
    numpy.savetxt('Outputs/' + geneName + '_NBTest.csv', pdat, fmt = ('%s'), delimiter = ',')
    numpy.savetxt('Outputs/' + geneName + '_NBFit.csv', fit, fmt = ('%s'), delimiter = ',')
    
def GenerateTransposonRecognitionSites(geneName, gene_fwd, gene_rev):
    print('Generating all possible five and eleven bp long fragments...')
    
    fragments_F = []
    fivemer_F = []
    elevenmer_F = []

    #fill forward frame
    for seq in [gene_fwd]:
        for i in range(0,len(seq)):
            if i < 34:
                index_adjustment = len(seq)-(34-i)
                sample = seq[index_adjustment:len(seq)] + seq[0:i+35]
                fragments_F.append(''.join(sample))
            if 34 <= i < (len(seq)-34):
                sample = seq[i-34:i+35]
                fragments_F.append(''.join(sample))
            if i >= (len(seq)-34):
                index_adjustment = 34-(len(seq)-i)+1
                sample = seq[i-34:len(seq)] + seq[0:index_adjustment]
                fragments_F.append(''.join(sample))
           
            #Pull out central 5mer and 11mers to match to corresponding transposon insertion frequencies 
            if i < 2:
                index_adjustment = len(seq) - (2-i)
                sample = seq[index_adjustment:len(seq)] + seq[0:i+2+1]
                fivemer_F.append(''.join(sample)) 
            if 2 <= i < (len(seq)-2):
                sample = seq[i-2:i+3]
                fivemer_F.append(''.join(sample)) 
            if i >= (len(seq)-2):
                index_adjustment = (2-(len(seq)-i))
                sample = seq[i-2:len(seq)] + seq[0:index_adjustment+1]
                fivemer_F.append(''.join(sample)) 
            if i < 5:
                index_adjustment = len(seq) - (5-i)
                sample = seq[index_adjustment:len(seq)] + seq[0:i+5+1]
                elevenmer_F.append(''.join(sample)) 
            if 5 <= i < (len(seq)-5):
                sample = seq[i-5:i+5+1]
                elevenmer_F.append(''.join(sample)) 
            if i >= (len(seq)-5):
                index_adjustment = (5-(len(seq)-i))
                sample = seq[i-5:len(seq)] + seq[0:index_adjustment+1]
                elevenmer_F.append(''.join(sample)) 


    #fill reverse frame
    fragments_NF = []
    fivemer_NF = []
    elevenmer_NF = []

    #fill forward frame
    for seq in [gene_rev]:
        for i in range(0,len(seq)):
            
            #69mer centered on residue 
            if i < 34:
                index_adjustment = len(seq)-(34-i)
                sample = seq[index_adjustment:len(seq)] + seq[0:i+35]
                fragments_NF.append(''.join(sample))
            if 34 <= i < (len(seq)-34):
                sample = seq[i-34:i+35]
                fragments_NF.append(''.join(sample))
            if i >= (len(seq)-34):
                index_adjustment = 34-(len(seq)-i)+1
                sample = seq[i-34:len(seq)] + seq[0:index_adjustment]
                fragments_NF.append(''.join(sample))
           
            #5mer cenetered on residue 
            if i < 2:
                index_adjustment = len(seq) - (2-i)
                sample = seq[index_adjustment:len(seq)] + seq[0:i+2+1]
                fivemer_NF.append(''.join(sample)) 
            if 2 <= i < (len(seq)-2):
                sample = seq[i-2:i+3]
                fivemer_NF.append(''.join(sample)) 
            if i >= (len(seq)-2):
                index_adjustment = (2-(len(seq)-i))
                sample = seq[i-2:len(seq)] + seq[0:index_adjustment+1]
                fivemer_NF.append(''.join(sample)) 
           
            #11mer centered on residue 
            if i < 5:
                index_adjustment = len(seq) - (5-i)
                sample = seq[index_adjustment:len(seq)] + seq[0:i+5+1]
                elevenmer_NF.append(''.join(sample)) 
            if 5 <= i < (len(seq)-5):
                sample = seq[i-5:i+5+1]
                elevenmer_NF.append(''.join(sample)) 
            if i >= (len(seq)-5):
                index_adjustment = (5-(len(seq)-i))
                sample = seq[i-5:len(seq)] + seq[0:index_adjustment+1]
                elevenmer_NF.append(''.join(sample)) 
    
    #reverse the reverse compliment
    fragments_NF = fragments_NF[::-1]
    elevenmer_NF = elevenmer_NF[::-1]
    fivemer_NF = fivemer_NF[::-1]          
    
    #Merge and adjust index to match VARIABLE.fasta, there is a more elegant fix for this
    fragments = fragments_F[5:len(fragments_F)] + fragments_F[0:5] + fragments_NF[5:len(fragments_F)] + fragments_NF[0:5]  
    fivemer =  fivemer_F[5:len(fivemer_F)] + fivemer_F[0:5] + fivemer_NF[5:len(fivemer_NF)] + fivemer_NF[0:5]
    elevenmer = elevenmer_F[5:len(elevenmer_F)] + elevenmer_F[0:5] + elevenmer_NF[5:len(elevenmer_NF)] + elevenmer_NF[0:5]

    print('\t'+ 'Number of Fragments (Fwd & Rev): ' + str(len(fragments)))
    print('\t'+ 'Number of Unique 5mers: ' + str(len(numpy.unique(fivemer))))
    print('\t'+ 'Number of Unique 11mers: '  + str(len(numpy.unique(elevenmer))))

    output_filename = 'Outputs/'+ geneName + '_PossibleSites.csv'
    rows = zip(range(1, len(fragments)+1), fragments, fivemer, elevenmer)
    with open(output_filename, 'wb') as f:
        writer = csv.writer(f)
        writer.writerow(['index', 'oligo', '5mer', '11mer'])
        for row in rows:
            writer.writerow(row)
            
    print('Written to CSV')
    print('\t' +  output_filename)
    print('\n')

def array(file):
    """
    Imports a csv of DNA sequences and inserts them into an array
    """
    sequences = []
    recSite = []
    freq = []
    with open(file, 'r') as csv_file:
        fileReader = csv.reader(csv_file, delimiter = "|")
        fileReader.next() # throwaway header row

        for row in fileReader:
            strippedRow = row[0].strip(",").split(',')
            sequences.append(strippedRow[1])
            recSite.append(strippedRow[2])
            freq.append(int(strippedRow[4]))

    return sequences, recSite, freq
