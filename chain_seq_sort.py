#!/usr/bin/env python3

from os import listdir
from os.path import isfile, join
from Bio.PDB.PDBParser import PDBParser
parser = PDBParser(PERMISSIVE=1, QUIET=True)


def filenames(files_dir):   # creates list object from directory
    names = [f for f in listdir(files_dir) if isfile(join(files_dir, f))]
    return names


with open('amino_acids.txt') as file:
    lines = file.readlines()
    acids_ls = [line.rstrip() for line in lines]  # creates single string in list from file
amino_acids = {}
for acid in acids_ls:
    acid = acid.split(',')
    amino_acids[acid[0].upper()] = acid[1]


def count_ch_sortRNA_RNP(file_path, dest_RNA, dest_RNP):
    #structure_id = '3BBI'
    #file_path = "path/to/3BBI.pdb"
    errors = []
    counter = 0
    count_RNA = 0
    count_RNP = 0
    RNP_mods = ['UNK', 'MSE', 'SEP', 'BB9', 'DHA', 'TPO', 'DBU', 'MH6', 'DCY', 'TS9', 'MLZ', 'ALY', 'MEA', 'SMC', 'HYP', 'ILX', 'TRX', 'CSX', 'CAS', 'DPR', 'DAB', '4J5', 'SCY', '0TD', 'OCS', 'D2T', 'MEQ', '2RX', 'CME', '4D4']
    prot_res = list(amino_acids.keys())
    prot_res.extend(RNP_mods)
    id_list = filenames(file_path)
    for id in id_list:
        if id[-4:] != '.pdb':
            continue
        structure_id = id[0:4]
        filename = file_path + '/'+ id

        try:
            structure = parser.get_structure(structure_id, filename)
            chain_seq = {}
            res_all = []
            count_protch = 0
            count_rnach = 0
            counter = counter + 1

            for chains in structure.get_chains(): #get chain_ID from structure
                chain_ID = str(chains).split('=')
                chain_ID = chain_ID[1].removesuffix('>')
                res_in_chain = []
                og_chain = []
                if chain_ID in chain_seq.keys(): #get just the first verion of sequence for given chain
                    continue
                else:
                    for res in chains.get_residues(): #get residues in sequence from chain_ID
                        res = str(res).split(' ')
                        #res_in_chain.append(res[1]) #appends residues for this chain_ID
                        res_all.append(res[1]) #appends to list of all residues together
                        if res[1] not in amino_acids.keys():
                            res_in_chain.append(res[1]) #appends residues for this chain_ID
                            og_chain.append(res[1])
                        else:
                            sngl_let_code = amino_acids[res[1]]
                            res_in_chain.append(sngl_let_code)
                            og_chain.append(res[1])

                chain_seq[chain_ID] = res_in_chain #add chain_ID:[sequence of res] to dict
                ch_idnty = next((acids for acids in prot_res if acids in og_chain),'NoAminoAcidRes')
                if ch_idnty == 'NoAminoAcidRes':
                    count_rnach = count_rnach + 1
                else:
                    count_protch = count_protch + 1

            identity = next((acids for acids in prot_res if acids in res_all),'NoAminoAcidRes') #check for first instance of amino acid in all residues
            if identity == 'NoAminoAcidRes':
                count_RNA = count_RNA+1
                with open(dest_RNA, 'a') as f:
                    f.write('#'+structure_id+'\t'+str(count_rnach)+'\t'+str(count_protch)+'\n')#adds line with dna/rna and protein chains
                    for key, value in chain_seq.items():
                        s = ".".join(map(str, value))
                        f.write(structure_id+'\t'+'%s\t%s' % (key, s)+'\n')#saves RNA dict as .txt (4P9R	A	C.A.U...)
                print('saved as RNA',structure_id,counter,'RNA/DNA:',count_rnach,'prot:', count_protch)
            else:
                count_RNP = count_RNP+1
                with open(dest_RNP, 'a') as f:
                    f.write('#'+structure_id+'\t'+str(count_rnach)+'\t'+str(count_protch)+'\n')#adds line with dna/rna and protein chains
                    for key, value in chain_seq.items():
                        v = ".".join(map(str, value))
                        f.write(structure_id+'\t'+'%s\t%s' % (key, v)+'\n')#saves RNP dict as .txt (4P9R	A	C.A.U...)
                print('saved as RNP',structure_id,counter,'RNA/DNA:',count_rnach,'prot:', count_protch)
        except UnicodeDecodeError:
            print('UnicodeDecodeError:',id)
            errors.append(id)
        except ValueError:
            print('ValueError', id)
            errors.append(id)
    print(counter,"RNA:",count_RNA, "RNP:", count_RNP)
    return errors

count_ch_sortRNA_RNP('file_path' , 'destination\RNA\.txt', 'destination\RNP\.txt')