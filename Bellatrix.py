#Current Time: 2 minutes for 20 proteins 
#Averages 6 seconds on each protein 

#Welcome to BELLATRIX
#Language: Python
#Note that you need the packages:
#numpy, urllib, tkinter, csv, time, os, pandas and biopandas to run Bellatrix. 
#Assuming you have pip installed, type pip install name_of_package into your terminal or command prompt where
#name_of_package is the name of the package you want to install. If pip isn't installed, please refer to online forums

#Basic User Instructions: 
#You have 2 options to select what protiens to create Star files of:
    #1 Input into the textbox and select 'enter'
    #2 Select a file: This file can either be a text file containing a list of proteins OR it can be a custom
    #.pdb file. Be sure to save it as a .pdb file before selecting the file! Note that you don't have to hit enter 
    
"""Once a file is picked, select 'begin construction'. This begins to create a star file library in the 
   local storage of your program. Read the output window or terminal (whatever your print() command displays to) to 
   see the progress of the creation of the library. Note that it will tell you if the file is missing atoms or residues."""

"""Once the library has been made, click the 'write to csv' to write it to a csv file. You will be prompted
   to create a name for this csv file, which will be saved to your working directory. Your CSV file will contain
   all of your star files. See https://2020.igem.org/Team:Calgary/Bellatrix for the anatomy of a Star file. Enjoy!"""
    

from biopandas.pdb import PandasPdb
import pandas as pd
import numpy as np 
import os
import time
import urllib
import tkinter as tk
import csv 


#this process uses the biopython package to download a pdb file, and then uses the
# biopandas package to download the PDB file

#PURPOSE: To get a list of the chains in the protein
#INPUT: A file like 'alphacarbons' or the original ATOMfile
#OUTPUT: A list of the chains sorted alphabetically
def getlistofchains (ATOMfile):
    chainlist = set(ATOMfile['chain_id'])
    chainlist = list(chainlist)
    chainlist = sorted(chainlist)
    print('The chain(s):')
    print(chainlist)
    return chainlist

    #this function checks if the residue numbers in the dataframe match the advertised number of residues
    #INPUT: original dataframe and alphacarbon dataframe
    #MODIFICATIONS: issues arise when there are skips in the data (missing atoms or residues)
    #OUTPUT: Returns True if numbers match and False otherwise
def checkresidues (alphacarbons, originalATOMfile):
    maxresidue = 0
    #first check and see if there are multiple chain ID's
    chainnumber = set(originalATOMfile['chain_id'])
    
    chainnumber = getlistofchains(originalATOMfile)
    #if there are multiple chain ID's, then we must find the max for each
    if (len(chainnumber)>1):
        print('There are multiple chains!')
        for i in range(len(chainnumber)):
            newmax = originalATOMfile.loc[originalATOMfile['chain_id'] == chainnumber[i]]
            newmax = newmax['residue_number'].max()
            maxresidue = maxresidue + newmax
    #otherwise, just get the max
    else:
        maxresidue = originalATOMfile['residue_number'].max()
        
    
    for i in range(len(chainnumber)):
        
        chainchunk = alphacarbons.loc[alphacarbons['chain_id']== chainnumber[i]]
        chainchunk.reset_index(inplace = True , drop = True )
        firstval = chainchunk.at[0, 'residue_number']
        maxresidue = maxresidue - firstval
    maxresidue = maxresidue + len(chainnumber)
        
    if (len(alphacarbons)== (maxresidue)):
       print('The Number of Residues Match')
       return True
    
    
    else:
        print('Check if there are other occupancies or missing residue numbers in file. The residue numbers dont match') 
        return False
    
#PURPOSE: When atoms positions are not 100% verified, there will be more than 1 entry for 
#the same atom. This causes the protein to appear to have too many residues. This function
#deletes all occupancies except the highest one
#INPUT: Dataframe containing the atomic coordinates
#OUTPUT: Dataframe with all alternate occupancies deleted except for the highest one
def easy_remove_low_occupancies (atoms):
    atoms = atoms.loc[(atoms['alt_loc'] !='B') & (atoms['alt_loc'] != 'C') & (atoms['alt_loc'] != 'D')]
    return atoms



#PURPOSE: To get a list of all subtraction matrices with distance vectors in an entire protein. The goal is to make
#this as efficient as possible
#INPUT: A pandas dataframe of alphacarbons
#OUTPUT:A subtraction matrix, where each row represents a specific alphacarbon's distance vectors between all other alphacarbons
# in the protein. It also returns a sequence dictionary that corresponds to the rows of the matrix. 
def get_all_of_submatrix(alphacarbons):
    alphacarbons.reset_index(inplace = True , drop = True )
    #this is for corresponding rows to residue numbers in matrix
    indexer = []
    dimensions = len(alphacarbons)
    matrix = np.zeros(shape = (dimensions, dimensions, 3))
    sequence = {}
    for i in range(dimensions): #row
        #coordinate at index i
        x = alphacarbons.at[i, 'x_coord']
        y = alphacarbons.at[i, 'y_coord']
        z = alphacarbons.at[i, 'z_coord']
        
        residue_name = alphacarbons.at[i, 'residue_name'] 
        residue_number = str(alphacarbons.at[i, 'residue_number']) + alphacarbons.at[i, 'chain_id']
        indexer.append(residue_number)
        sequence[residue_number] = residue_name
        for j in range(dimensions): #column
            #coordinate at index j
            if (j<i):
                matrix[i][j]= matrix[j][i] * -1 
                continue
            x1 = alphacarbons.at[j, 'x_coord']
            y1 = alphacarbons.at[j, 'y_coord']
            z1 = alphacarbons.at[j, 'z_coord']
            if (i == j):
                matrix[i][j]= (0, 0, 0)

            else:
                x2 = x1-x
                y2 = y1-y
                z2 = z1-z
                matrix[i][j][0]=x2
                matrix[i][j][1]=y2
                matrix[i][j][2]=z2
                
    return matrix, sequence, indexer

#this function will determine what experimental methods were used to get the pdb file. 
#INPUT: file path to the PDB
#OUTPUT: string of the experimental data type
def findEXPtype(pdbfilepath):
    fp = open(pdbfilepath, "r")
    experiment_type = 'NODATA'
    while (1):
        c = fp.readline()
        c = c.split()
        if (c[0] == 'EXPDTA'):
            experiment_type = ''
            for i in range(1, len(c)):
                experiment_type = experiment_type + c[i]
            fp.close()
            return experiment_type
        #this is in the event that there is no input for EXPDTA
        elif (c[0]=='AUTHOR'):
            fp.close()
            break

#PURPOSE: Read the header of the PDB file and return a string containing a warning of the types of missing data 
#NOTES:
#Be aware that REMARK 470 in a PDB file is were it lists the missing atoms 
def check_for_missing_data(pdbfilepath):
    fp = open(pdbfilepath, "r")
    errortype = ""
    while (1):
        c = fp.readline()
        if('DBREF' in c):
            fp.close()

            return errortype

        if('MISSING' in c):
            pos = c.find('MISSING')
            errortype = errortype + " " + c[pos:len(c)] + '\0'

#PURPOSE: To get the missing residues so as to display this in our library as part of the metadata
#INPUT: Filepath to the pdb file being analyzed 
#OUTPUT: a list of the missing residues
#MODIFICATIONS: This process can be expedited  and streamline process by knowning what is missing
def get_missing_residues(pdbfilepath):
    fp = open(pdbfilepath, "r")
    missing_residues = []
    while(1):
        c = fp.readline()
        if ('REMARK 465' in c): 
            while(1):
                c = fp.readline()
                if (not 'REMARK 465' in c):
                    fp.close()
                    return missing_residues
                
                if ('M RES C SSSEQI' in c):
                    while(1):
                        c = fp.readline()
                        if (not 'REMARK 465' in c):
                            break
                        c = c.split()
                        missing_residues.append(c[2:len(c)])
        elif('DBREF' in c):
            fp.close()
            return missing_residues

        
#tkinter stuff ------------------------------------------------------

def read_in_text():
    global cellulases
    cellulases = e1.get().split(',')
    
#As long as the input is separated by white space, you can do it
def read_in_text_file():
    from tkinter.filedialog import askopenfilename
    text_file_filepath = askopenfilename()
    global cellulases
    #if the file the user chooses is of .pdb type, then get rid of cellulases and just read that one file. The first 
    #entry of cellulases is the CUSTOM indicator, and the second entry is the filepath to the .pdb file
    if '.pdb' in text_file_filepath:
        cellulases = ['CUSTOM', text_file_filepath]
        return 
    
    
    with open(text_file_filepath, "r") as file:
        holder = []
        list_from_text_file = []
        c = file.readlines()
        for x in c: 
            x = x.split()
            if (x != []):
                holder.append(x)
        for x in holder:
            for i in range (len(x)):
                list_from_text_file.append(x[i])

        cellulases = list_from_text_file
    
def begin_bellatrix (): 
    
    
    global cellulases, mylib, NMR, nonexistent
    

    #library that you are storing star files in
    mylib ={} 
    
    start = time.time()
    if cellulases[0] != 'CUSTOM':
        number_of_proteins =  len(cellulases)

        #list of structures that are nonexistent 
        nonexistent = []
        NMR = []
        
        for i in range (number_of_proteins):

            print("Current Progress: %(i)s/%(total)s" %{'i': i+1, 'total':number_of_proteins})
            #the protein PDB code lol
            pteen = cellulases[i]

            #this attempts to download the pdb file from the Protein Data Bank
            print('Beginning download for ' + pteen)
            try: 
                urllib.request.urlretrieve('https://files.rcsb.org/download/' + pteen + '.pdb', pteen + '.pdb')

            except urllib.error.HTTPError:
                print('The pdb file of ' + pteen + ' could not be downloaded. It might not exist.')
                print('')


            #this creates the biopandas structure that easily converts PDB files into a Pandas dataframe
            ppdb = PandasPdb()


            #gets the working directory to read the file
            urdirectory = os.getcwd()

            #this is the path to our downloaded file
            filepath = '%(currentdirec)s/%(code)s.pdb' %{'currentdirec': urdirectory, 'code': pteen}




            #read the pdb file into the biopandas dataframe, but also use exception
            try: 
                ppdb.read_pdb(filepath)

            except FileNotFoundError: 
                #this records the nonexistent file
                nonexistent.append(pteen)
                continue

            except: 
                print('Other error in opening the file of' + pteen)


            exptype = findEXPtype (filepath)
            print(exptype)

            #handle NMR files
            if('NMR' in exptype):
                print("The protein " + pteen +" could not be added because it is an NMR file")
                NMR.append(pteen)
                os.remove(filepath)
                continue

            missing_info = check_for_missing_data(filepath)
            print(missing_info)
            #this gets the ATOM portion of the PDB and loads it into the dataframe. 
            atomstuff = ppdb.df['ATOM']

            #this gets the alphacarbons into one dataframe
            alphacarbons = atomstuff.loc[atomstuff['atom_name'] == 'CA']


            #this portion handles the errors one might find in the PDB file. It always parses the file correctly (to date no issues), 
            #and if it says it doesnt, it means that the file has some issues in itself.(Missing residues or atoms)
            check = checkresidues(alphacarbons, atomstuff)
            if (check == False):
                alphacarbons = easy_remove_low_occupancies(alphacarbons)
                check = checkresidues(alphacarbons, atomstuff)
            if (check == False):
                print('Occupancies are not the issue, check the file for errors for %s. (Could be skips in residue numbers)' %pteen)

            matrix, sequence, indexer = get_all_of_submatrix(alphacarbons)

            missing_residues = get_missing_residues(filepath)
            print('The number of missing residues: %s' %len(missing_residues))


            
            #create a dictionary of entries, allowing the user to manipulate
            # the data in their program memory, or they canformat it into a CSV file.


            mylib[pteen] = {'stars': [matrix, indexer], 'residues': sequence, 'missing': missing_residues, 'metadata':[missing_info, exptype]}

            print()
            #this deletes the pdb file that you were working with 
            os.remove(filepath)
            
       


        #now modify the list of proteins so that you print the correct list to the csv 
        for i in range (len(nonexistent)):
            cellulases.remove(nonexistent[i])

        for i in range (len(NMR)): 
            cellulases.remove(NMR[i])

    else:
            #this is what you do if cellulases = "custom" (you are reading in a custom file)

        print('Beginning of Bellatrix Creation on Custom File: ' + cellulases[1])
        ppdb = PandasPdb()
        filepath = cellulases[1]
        ppdb.read_pdb(filepath)
        atomstuff = ppdb.df['ATOM']
        #this checks for chain ID's
        alphacarbons = atomstuff.loc[atomstuff['atom_name'] == 'CA']
        if atomstuff.at[0, 'chain_id'] == '' or atomstuff.at[5, 'chain_id'] == '':
            start_of_chain = alphacarbons.loc[alphacarbons['residue_number']==1] 

            #this contains the index of the start of each chain 
            start_of_chain = list(start_of_chain.index.values)

            #now we need to put in the chains into each row in the data frame 
            #get indexes of the entire pandas dataframe
            main_indexer = list(alphacarbons.index.values)
            #this loop goes through the alphacarbons and inputs the correct chain ID's
            chain_input = 'A'
            main_indexer_i = 0
            for i in main_indexer:
                #if the index equals the start of a new chain, change what you set the chain ID to
                if main_indexer_i<len(start_of_chain) and i == start_of_chain[main_indexer_i]:
                    chain_input = chr(ord(chain_input) + 1)
                    main_indexer_i = main_indexer_i + 1

                alphacarbons.at[i, 'chain_id'] = chain_input

        #if there are no issues with chain ID
        print('starting to create star file')
        matrix, sequence, indexer = get_all_of_submatrix(alphacarbons)
        print('Star file successfully created')
        missing_residues = 'NO MISSING RESIDUES'
        missing_info = 'NO MISSING INFO, CUSTOM FILE'
        exptype = 'Experimental Type: Dynamic Data'
        mylib[1] = {'stars': [matrix, indexer], 'residues': sequence, 'missing': missing_residues, 'metadata':[missing_info, exptype]}

    finish = time.time()
    print(f"\n Bellatrix Library created in {finish - start} seconds")


            
    
#PURPOSE: To write the Bellatrix library to a csv file once the button on the GUI is pressed
#INPUT: Nothing, but this is linked to the button on the GUI
#OUTPUT: Nothing, writes Bellatrix library to Star file
def write2csv():
    #create a list of dictionaries for every row of csv file
    print('Beginning to write library to a CSV file ')
    start = time.time()
    gucci_list = []
    counter = 0
    
    #this section gets user input for csv file name
    from tkinter import simpledialog

    ROOT = tk.Tk()

    ROOT.withdraw()
    # the input dialog
    USER_INP = simpledialog.askstring(title="Bellatrix",
                                  prompt="Please input what you would like to name your csv file (without .csv):")
    csv_name = USER_INP + '.csv'
    with open(csv_name, 'w', newline = '') as file:
        #ADD IN CODE TO PRINT OUT LIBRARY FIRST 

        writer2 = csv.writer(file)
        writer2.writerow(cellulases)

        for key in mylib:
            counter = counter + 1
            length = len(mylib[key]['stars'][1])
            print("Current Progress: %(i)s/%(total)s" %{'i': counter, 'total':len(cellulases)})
            #get the metadata and protien name into a list to write to CSV file
            firststuff = []
            firststuff.append(key)
            firststuff.append(mylib[key]['metadata'])
            firststuff.append(mylib[key]['missing'])

            for i in range(length):
                numbers = {}
                numbers['residue_number'] = mylib[key]['stars'][1][i]
                numbers['amino_acid'] = mylib[key]['residues'][mylib[key]['stars'][1][i]]
                fieldnames =['residue_number', 'amino_acid']
                for j in range(length):
                    numbers[j] = mylib[key]['stars'][0][i][j] #might have to make this as a tuple? 
                    fieldnames.append(j) #try and simplify this and take the last persons code



                writer = csv.DictWriter(file, fieldnames = fieldnames)
                writer2 = csv.writer(file)
                if(i == 0):
                    writer2.writerow(firststuff)
                    writer.writeheader()
                writer.writerow(numbers)
            
            print(key +' has been converted to a csv format')
    end = time.time()
    print(f"Runtime of the program is {end - start}")
    

#------------------------------------------ Body of Code that runs everything  --------------------------------



myprog = tk.Tk()
myprog.title('BELLATRIX')

#variables in use
cellulases = ""
mylib = {}
var1 = tk.IntVar(myprog)
chi_angle_check = 0



#get the inputted protein list
#This reads in any text you have inputted into the GUI.

tk.Label(myprog, text = 'Enter your desired list of proteins separated by a comma').grid(row=0, column = 0)
user_input = tk.StringVar(myprog)
#now add a box for entries
e1 = tk.Entry(myprog)
e1.grid(row = 0, column=1)

#get a button to allow the program to read in the text
enter_button = tk.Button(myprog, text = 'Enter', width = 20, command =read_in_text)
enter_button.grid(row = 0, column = 2)

#button to offer the option to upload a list of proteins via text file 
select_file_button = tk.Button(myprog, text = 'Select File', width = 20, command =read_in_text_file)
select_file_button.grid(row = 0, column = 3)



#get a button to allow the program to read in the text
enter_button = tk.Button(myprog, text = 'Enter', width = 20, command =read_in_text)
enter_button.grid(row = 0, column = 2)

#button to offer the option to upload a list of proteins via text file 

select_file_button = tk.Button(myprog, text = 'Select File', width = 20, command =read_in_text_file)
select_file_button.grid(row = 0, column = 3)

#button to begin Bellatrix
button = tk.Button(myprog, text = 'Begin Construction', width = 25, command =begin_bellatrix)
button.grid(row = 5, column = 0)

#button to write to csv
csv_button = tk.Button(myprog, text = 'Write to CSV file', width = 25, command =write2csv)
csv_button.grid(row = 6, column = 0)




myprog.mainloop()



