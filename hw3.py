def alignment(input_path, score_path, output_path, aln, gap):
    #LOAD PACKAGES
    import numpy as np
    import pandas as pd
    import Bio as bio

    #READ DATA
    score_mat = pd.read_csv(score_path, header = 0, comment = "#", sep = "\s+")
    seq = pd.read_csv(input_path, header = None)
    #Extract protein sequences as seq1 and seq2
    seq1 = pd.DataFrame(list(seq.at[1,0]))
    seq2 = pd.DataFrame(list(seq.at[3,0]))
    #Extract sequence ids (the names of the proteins) as protein1 and protein2
    protein1 = str(seq.at[0,0])
    protein1 = protein1[1:]
    protein2 = str(seq.at[2,0])
    protein2 = protein2[1:]

    #SEQUENCE ALIGNMENT
    #Add an extra row & extra column to the alignment table (alg_table) in case we start with a gap.
    alg_table = np.zeros(shape = (len(seq1) + 1, len(seq2) + 1))
    #Fill the top row & leftmost column.
    if aln == "global":
        alg_table[:,0] = np.arange(len(seq1) + 1) * gap
        alg_table[0,:] = np.arange(len(seq2) + 1) * gap
    else: 
        alg_table[:,0] = np.zeros(len(seq1) + 1)
        alg_table[0,:] = np.zeros(len(seq2) + 1)

    #Create a table to keep track of the direction the alignment is moving

    #0: the alignment won't move to this cell
    #1: the alignment came from the top left cell (a match/mismatch)
    #2: the alignment came from the left cell
    #3: the alignment came from the top cell

    d_table = np.zeros(shape = (len(seq1) + 1, len(seq2) + 1))
    if aln == "global":
        d_table[:,0] = np.ones(len(seq1) + 1) * 3
        d_table[0,:] = np.ones(len(seq2) + 1) * 2
        d_table[0,0] = 0
    else:
        d_table[:,0] = np.zeros(len(seq1) + 1) 
        d_table[0,:] = np.zeros(len(seq2) + 1)
    
    #Fill in the tables
    for i in range(1, len(seq1) + 1):
        for j in range(1, len(seq2) + 1):

            match = alg_table[i-1,j-1] + score_mat.at[seq1.at[i-1,0], seq2.at[j-1,0]]
            right = alg_table[i,j-1] + gap
            down = alg_table[i-1,j] + gap
            
            if aln == "global":
                chosen = max(match, right, down)
            else: chosen = max(match, right, down, 0)
            
            alg_table[i, j] = chosen
            if chosen == match:
                d_table[i, j] = 1
            elif chosen == right:
                d_table[i, j] = 2
            elif chosen == down:
                d_table[i, j] = 3
            else:
                d_table[i, j] = 0
    #TRACEBACK
    if aln == "local":
        #For local alignment, find the locations with the highest score on alg_table
        max_position = np.where(alg_table == np.max(alg_table, axis=None))
        #The traceback(s) corresponding to the highest scores 
        traceback_list = []
        

        for i in range(len(max_position[0])):
            latest_row = max_position[0][i]
            latest_column = max_position[1][i]
            previous_action = d_table[latest_row, latest_column]
            traceback = [previous_action]
            while previous_action != 0:
                if previous_action == 1:
                    latest_row = latest_row - 1
                    latest_column = latest_column - 1
                elif previous_action == 2:
                    latest_column = latest_column - 1
                else:
                    latest_row = latest_row -1
                previous_action = d_table[latest_row, latest_column]
                traceback.append(previous_action)
            traceback.reverse()
            traceback.pop(0)
            traceback_list.append(traceback)
            
        #The traceback(s) with the highest score and longest length are selected
        traceback_selected = []
        for i in range(len(traceback_list)):
            if len(traceback_list[i]) == len(max(traceback_list, key=len)):
                traceback_selected.append(traceback_list[i])
    else:
        #For global alignments
        previous_action = d_table[len(seq1), len(seq2)]
        latest_row = len(seq1)
        latest_column = len(seq2)
        traceback = [previous_action]

        while previous_action != 0:
            if previous_action == 1:
                latest_row = latest_row - 1
                latest_column = latest_column - 1
            elif previous_action == 2:
                latest_column = latest_column - 1
            else:
                latest_row = latest_row -1
            previous_action = d_table[latest_row, latest_column]
            traceback.append(previous_action)
        traceback.reverse()
        traceback.pop(0)
        traceback_selected = [traceback]
    
    #CREATE ALIGNED SEQUENCES
    aseq1_list = []
    aseq2_list = []

    for t in range(len(traceback_selected)):
        traceback = traceback_selected[t]
        aseq1 = []
        aseq2 = []
        aseq1_index = latest_row
        aseq2_index = latest_column
        for i in range(len(traceback)):
            if traceback[i] == 1:
                aseq1.append(seq1.at[aseq1_index,0])
                aseq2.append(seq2.at[aseq2_index,0])
                aseq1_index = aseq1_index + 1
                aseq2_index = aseq2_index + 1
            elif traceback[i] == 2:
                aseq1.append("-")
                aseq2.append(seq2.at[aseq2_index,0])
                aseq2_index = aseq2_index + 1
            else: 
                aseq1.append(seq1.at[aseq1_index,0])
                aseq2.append("-")
                aseq1_index = aseq1_index + 1
        aseq1_list.append(aseq1)
        aseq2_list.append(aseq2)
    
    #WRITE FASTA FILE
    #Source for following code: https://milliams.com/courses/biopython/Input%20and%20Output.html

    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio import SeqIO

    records = []

    for i in range(len(aseq1_list)):
        rec1 = SeqRecord(
            Seq(''.join(aseq1_list[0])),
            id = protein1,
            description="")
        rec2 = SeqRecord(
            Seq(''.join(aseq2_list[0])),
            id = protein2,
            description="")
        records.append(rec1)
        records.append(rec2)

    SeqIO.write(records, output_path, "fasta")



