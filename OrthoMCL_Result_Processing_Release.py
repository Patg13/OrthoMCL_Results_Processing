#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# DEPENDENCIES
# matplotlib 3.3.4 (Required)
# Spacy 3.0.3 (Optionnal, use --ignore_spacy flag to skip)
# en_core_web_lg dictionnary for Spacy

import sys,os
import argparse
import itertools
import matplotlib
import re
matplotlib.use("Agg") # Non interactive Backend of Matplotlib


# venn_l is a python3 script that's on the HDD with this script
# I am not the creator of this script, original reference: https://github.com/tctianchi/pyvenn
import venn_l as venn

print("OrthoMCL Result formatting Program V0.2 ALPHA REL")
print("By: Patrick Gagne")

# TODO Use better comment to define objects (look for correct format)

parser=argparse.ArgumentParser(description='Process OrthoMCL raw results into human interpretable data V0.2 ALPHA RELEASE')

parser.add_argument("-g","--group", dest="group_file",required=True, help=("Group file containing proteins groups [REQUIRED]"))
parser.add_argument("-s","--singleton", dest="singletons_file",required=False, help=("Singleton file containing single proteins  [OPTIONAL]"))
parser.add_argument("-p","--proteome", dest="proteome_file",required=False, help=("Fasta file regrouping all proteomes (MUST INCLUDE PROTEIN ANNOTATION AND ID MUST BE SEPERATE BY A SPACE) [OPTIONAL]"))
parser.add_argument("-c","--species", dest="spc_list", help=("Species tag list (ex: TBE,COVI,...)  [REQUIRED]"))
parser.add_argument("-r","--regroup", dest="regroup_list",required=False, help=("List of species you want to regroup (; to link ,to seperate) SP1;SP2,SP3;SP4 means 1 regroup with SP1 and SP2 and 1 regroup with SP3 and SP4  [OPTIONAL]"))
parser.add_argument("-l","--limit", dest="table_limit",type=int,default=11, help=("Ortho Table limit for cluster count (every cluster >= this value will have their count add to the same value (Val+)) ! higher value means more columns in the table ! [DEFAULT=11]"))
parser.add_argument("--venn_dpi", dest="dpi",type=int,default=400, help=("Pixel Density for Venn Diagram greater value means sharper but bigger png [DEFAULT=400]"))
parser.add_argument('--ignore_spacy', action='store_true', help=("Ignore Spacy step to reduce annotation count (remove spacy dependency, MUCH faster, but will produce a massive annotation CSV)"))
parser.add_argument("--spacy_simthres", dest="spacy_thres",type=float,default=0.95, help=("Similarity threshold for spacy annotation regrouping MODIFY AT YOUR OWN RISK [DEFAULT=0.95]"))
parser.add_argument('--spacy_pipeline', dest="spacy_pipe",default="en_core_web_lg", help=("Spacy language pipeline to use (all the tests were done with en_core_web_lg, so it's not recommended to use another) [default=en_core_web_lg]"))
parser.add_argument("-o","--output", dest="output_folder",default="orthomcl_resgen_out/", help=("Output folder name [Default=orthomcl_resgen_out/]"))

args=parser.parse_args()

if args.ignore_spacy == False and args.proteome_file != None:
    if args.spacy_thres > 1.0:
        print("ERROR: Spacy Similarity Threshold cannot be more then 1.0")
        sys.exit(1)
    print("Importing Spacy library + model")
    try:
        import spacy # Spacy is a text comparison library using statistical models
    except ModuleNotFoundError:
        print("ERROR: Spacy module not found, you must install it : pip install spacy")
        sys.exit(1)
    try:
        nlp = spacy.load(args.spacy_pipe) # en_core_web_lg MUST BE downloaded as a package with this command : py -m spacy download en_core_web_lg
    except OSError:
        print("ERROR: %s pipeline not found by spacy"%(args.spacy_pipe))
        print("if you didn't specified another pipeline, you must download the en_core_web_lg by using this command: py -m spacy download en_core_web_lg")
        sys.exit(1)
    
if args.output_folder[-1] != "/":
    args.output_folder+="/"
os.makedirs(os.path.dirname(args.output_folder), exist_ok=True)

# Function to remove illegal character in s (so s can be a filename)
def get_valid_filename(s):
    s = str(s).strip().replace(' ', '_')
    return re.sub(r'(?u)[^-\w.]', '', s)

# Create an ortho group that will contain the name of the group and the list of protein
# Will also contain certain stats that will be useful for other functions
class ortho_group:
    def __init__(self,name,protein_list):
        self.name=name
        self.species={}
        self.proteins={}
        # Create a dictionnary where the species name will be used as key and will contain all associated proteins
        # Will also keep a trace of which species has been found in the group and it's count
        for i in protein_list:
            spl=i.split("|")
            try:
                self.proteins[spl[0]].append(spl[1].replace("\n",""))
            except KeyError:
                self.proteins[spl[0]]=[spl[1].replace("\n","")]
            try:
                self.species[spl[0]]+=1
            except KeyError:
                self.species[spl[0]]=1
        self.species_count=len(self.species.keys())
        self.protein_count=len(protein_list)
        self.species_list=sorted(list(self.species.keys()))
        self.classification=self.classify_group()
    # Just in case we need to have a count of each species in the group in a list of tuple format instead of a dictionnary
    # The returned information is basically the same as self.species
    def count_by_species(self):
        result_list=[]
        for i in self.species_list:
            result_list.append((i,len(self.proteins[i])))
        return result_list
    # Each group fit in a category according to the ortho table, so this function is useful to determine that category
    def classify_group(self):
        if self.species_count == 1:
            # I use the name "branch" because the folder structure where this info is used is a tree
            master_branch="unique"
        else:
            master_branch="shared"
        if self.protein_count >= args.table_limit:
            sub_branch=master_branch+"_"+str(args.table_limit)+"+"
        else:
            sub_branch=master_branch+"_"+str(self.protein_count)
        return (master_branch,sub_branch)
    

# Will generate a complete count table to determine all the species combination
class ortho_table:
    def __init__(self,comb_dicti,limit):
        self.lineheader=comb_dicti.keys()
        self.proteins_count=[]
        for i in self.lineheader:
            self.proteins_count.append([i,sum(comb_dicti[i])])
        self.colheader=["Species"]
        count_total_dict={}
        for i in range(1,limit):
            self.colheader.append(str(i))
            count_total_dict[i]=0
        count_total_dict[limit]=0
        count_total_dict[limit+1]=0
        self.colheader.append(str(limit)+"+")
        self.colheader.append("Total")
        self.lines=[]
        
        for i in self.lineheader:
            #print(i)
            line=[i]
            #print(line)
            for j in range(1,limit):
                count=comb_dicti[i].count(j)
                line.append(count)
                count_total_dict[j]+=count
            total_len=len(comb_dicti[i])
            count_total_dict[limit]+=total_len-sum(line[1:])
            line.append(total_len-sum(line[1:]))
            count_total_dict[limit+1]+=sum(line[1:])
            line.append(sum(line[1:]))
            
            self.lines.append(line)
        # Under Total line calculation
        total_line=["Total"]
        for i in sorted(count_total_dict.keys()):
            total_line.append(count_total_dict[i])
        self.lines.append(total_line)
    # Directly print the table to the output (with basic formatting)
    def print_table(self,sep="\t"):
        print(sep.join(self.colheader))
        for i in self.lines:
            print(sep.join([str(elem) for elem in i]))
            
    def save_table_to_file(self,file,sep=","):
        # CSV format will be used by default (sep=",")
        file.write(sep.join(self.colheader)+"\n")
        for i in self.lines:
            file.write(sep.join([str(elem) for elem in i])+"\n")

    def prep_venn(self,mode):
        # This function will return a label dictionnary (no need to get label with original data, which would be much harder to do)
        # The size is determined automatically according to the number of lines in the table
        size_verif=[3,7,15,31,63] # Number of line necessary for 2,3,4,5,6 level venn diagram (7+ not supported)
        data=[]
        if len(self.lines[:-1]) not in size_verif:
            print("WARNING: Cannot generate venn diagram: inconsistent number of lines in ortho_table")
            print("This is normal if you have more than 6 species without regrouping (max level: 6)")
            return ""
        size=size_verif.index(len(self.lines[:-1]))+2
        print("Venn size detected : %d"%(size))
        # the library venn_l support venn = 6, but use triangle instead of ellipses, which is harder to read
        if size == 6:
            print("WARNING: size 6 venn diagram detected, this venn will be very hard to read")
            print("You should REALLY regroup at least 2 species together")
        # MODE Cluster (will use cluster count to generate the labels)
        if mode == "C":
            for i in self.lines[:-1]:
                data.append((i[0],i[-1]))
        # MODE Proteins (will use protein count to generate the labels)
        if mode == "P":
            for i in self.proteins_count:
                data.append((i[0],i[1]))
        # MODE Both (will use cluster count AND protein count to generate the labels)
        if mode == "B":
            for i in range(0,len(self.proteins_count)):
                data.append((self.lines[i][0],self.lines[i][-1],self.proteins_count[i][1]))

        if mode != "P" and mode != "C" and mode != "B":
            print("WARNING: Unrecognized venn diagram mode...skipping")
            return ""
        s_grp={}
        l_grp=[]
        count=0
        for i in data:
            if len(i[0].split("|")) == 1:
                l_grp.append(i[0])
                s_grp[i[0]]=count
                count+=1
        
        label_dict={}
        # Generate a binary list which will determine the identity of the venn diagramm groups
        # Then associate each count with it's own binary to generate the final dictionnary
        # The binary is the dict key and the count is the value
        for i in data:
            main_bin=["0"] * size
            for j in i[0].split("|"):
                main_bin[s_grp[j]]="1"
            if mode != "B":
                label_dict["".join(main_bin)]=i[1]
            else:
                # Cluster on top and protein count in parenthesis under
                label_dict["".join(main_bin)]=str(i[1])+"\n"+"("+str(i[2])+")"
        return (l_grp,label_dict)

        
species_list=sorted(args.spc_list.split(","))

# Species Regrouping Section (Will not be executed if no regrouping list is given)
if args.regroup_list != None:
    rep_dict={}
    print("DEBUG: Species Regrouping asked")
    print("DEBUG: Isolation lists")
    regroups=args.regroup_list.split(",")
    print("DEBUG: Verifying Regroup Validity")

    # We must test for the regroup validity
    test=args.regroup_list.replace(";",",").split(",")
    testset=set(test) # Use a set because it can't contain duplicate
    # If there is a disparity between the length of the list and the set, it means there is a duplicate in the regroup
    if len(test) != len(testset):
        print("ERROR: All species in regroup list can only be used once")
        sys.exit(2)
    # Here is to test for species not in the species list (typo for exemple)
    for i in test:
        if i not in species_list:
            print("ERROR: %s from regroup list not in your species list"%(i))
            sys.exit(2)
    # If there is a regroup with only one species
    for i in regroups:
        if len(i.split(";")) == 1:
            print("ERROR: Cannot regroup a single species")
            sys.exit(2)
    print("DEBUG: Regroup list valid, processing...")
    #Create a dictionnary where the key is the original name and the value is the new name
    for i in regroups:
        spl=i.split(";")
        for j in sorted(spl):
            rep_dict[j]=i.replace(";","+") # Change + here to change the seperator between regroups EX: SP1+SP2

    # We must also add the original species (those who are not regrouped) in the dictionnary
    # It will make the replacing easier
    for i in species_list:
        try:
            test=rep_dict[i]
        except KeyError:
            rep_dict[i]=i
    print("DEBUG: Adjusting Species list with groups")
    # We must also modify the species list to reflect the new groups (we use a set to remove duplicate)
    new_species=set()
    for i in species_list:
        try:
            new_species.add(rep_dict[i])
        except KeyError:
            new_species.add(i)
    species_list=sorted(list(new_species))
    

    


ortho_grps=[]
print("Processing Groups into ortho groups") 
with open(args.group_file,"r") as infile:
    for line in infile:
        spl=line.split(" ",1)
        if args.regroup_list == None:
            ortho_grps.append(ortho_group(spl[0],spl[1].replace("\n","").split(" ")))
        else:
            # These instructions will replace the species name for the regrouped ones
            prot_line=spl[1].replace("\n","").split(" ")
            new_line=[]
            for i in prot_line:
                sp=i.split("|")[0]
                # OLD LINE : new_line.append(i.replace(sp,rep_dict[sp]))
                # 1 was added in case the replacement string was also in the sequence name
                new_line.append(i.replace(sp,rep_dict[sp],1))
            ortho_grps.append(ortho_group(spl[0],new_line))


if args.singletons_file != None:
    print("Processing Singletons file into single ortho groups")
    s_count=0
    # Create an objet ortho_grp for each singletons the same way and give it an unique name (single_#)
    with open(args.singletons_file,"r") as infile:
        for line in infile:
            newname="single_"+str(s_count)
            if args.regroup_list == None:
                ortho_grps.append(ortho_group(newname,line.replace("\n","").split(" ")))
            else:
                # Same thing is done in groups files, it's just easier here because there is no protein list, only single ones
                sp=line.split("|")[0]
                ortho_grps.append(ortho_group(newname,line.replace("\n","").replace(sp,rep_dict[sp]).split(" "))) # Split must be done anyway to imitate the format ortho_grp needs
            s_count+=1

        
# Combination table (will generate every ortho groups possibilities) 
print("Processing Ortho groups into an ortho table")
comb_list=[]
for i in range(1,len(species_list)+1):
    for comb in itertools.combinations(species_list,i):
        comb_list.append(comb)


comb_dict={}
comb_list_sort=sorted(comb_list,key=lambda x:len(x),reverse=False)
for i in comb_list_sort:
    comb_dict["|".join(i)]=[]

for i in ortho_grps:
    comb_dict["|".join(i.species_list)].append(i.protein_count)

ortho_tb=ortho_table(comb_dict,args.table_limit)

# Write OrthoTable csv
otable_file=open(args.output_folder+"OrthoTable.csv",'w')
ortho_tb.save_table_to_file(otable_file,",")
otable_file.close()

# Producing the Protein count, the Cluster count and the combined Venn diagrams
print("Generating Venn diagrams")
for i in ["P","C","B"]:
    venn_data=ortho_tb.prep_venn(i)
    if venn_data !="":
        if len(venn_data[0]) == 2:
            fig, ax = venn.venn2(venn_data[1],names=venn_data[0])
            fig.savefig(args.output_folder+i+"_venn_L2.png",dpi=args.dpi)
            fig.clf()
            
        if len(venn_data[0]) == 3:
            fig, ax = venn.venn3(venn_data[1],names=venn_data[0])
            fig.savefig(args.output_folder+i+"_venn_L3.png",dpi=args.dpi)
            fig.clf()
           
        if len(venn_data[0]) == 4:
            fig, ax = venn.venn4(venn_data[1],names=venn_data[0])
            fig.savefig(args.output_folder+i+"_venn_L4.png",dpi=args.dpi)
            fig.clf()
            
        if len(venn_data[0]) == 5:
            fig, ax = venn.venn5(venn_data[1],names=venn_data[0])
            fig.savefig(args.output_folder+i+"_venn_L5.png",dpi=args.dpi)
            fig.clf()
            
        if len(venn_data[0]) == 6:
            fig, ax = venn.venn6(venn_data[1],names=venn_data[0])
            fig.savefig(args.output_folder+i+"_venn_L6.png",dpi=args.dpi)
            fig.clf()
            
    else:
        print("WARNING: Empty data to generate venn_diagram...skipping")



# NOW READY TO GENERATE FASTA FOLDERS (BY GROUPS AND ACCORDING TO ORTHO TABLE)

if args.proteome_file == None:
    print("WARNING: No proteome file has been provided...exiting program")
    sys.exit(0)

# This step will read the proteome file to associate ortho groups proteins to their annotations and amino acid sequences
print("Extracting annotation information and protein sequences from proteomes file")
anno_dict={}
seq_dict={}
warning=0
with open(args.proteome_file,"r") as infile:
    for line in infile:
        if line[0] != ">":
            # The use of a list as value here is just in case the proteome file is not in the single line fasta format
            try:
                seq_dict[spl[0].replace(">","")].append(line.replace("\n",""))
            except KeyError:
                seq_dict[spl[0].replace(">","")]=[line.replace("\n","")]
        else:
            spl=line.split(" ",1)
            # in case a protein in proteome file does not have an annotation, \n must be replaced here
            if len(spl) == 1:
                spl[0]=spl[0].replace("\n","")
            # NEW TRY in case a protein in proteome file does not have an annotation
            try:
                anno_dict[spl[0].replace(">","")]=spl[1].replace("\n","")
            except IndexError:
                warning=1
                anno_dict[spl[0].replace(">","").replace("\n","")]="NA"

if warning == 1:
    print("WARNING: Some sequences in the proteomes file does not have annotation")
                


# Creating a dictionnary for each file to write
# Using this format : SPECIES_LIST;X_MEMBER : [ list of ortho_groups ]
# gbfs = Genes By Family By Species
print("Generating Genes by families by species csv files")
gbf_dict={}
# Classification of each ortho group depending on its species list and its protein count 
for i in ortho_grps:
    # Skipping Singletons, 1 species list and less than 3 proteins
    if len(i.species_list) == 1:
        continue
    if i.protein_count == 2:
        continue
    try:
        gbf_dict["_".join(i.species_list)+";"+str(i.protein_count)+"_members"].append(i)
    except KeyError:
        gbf_dict["_".join(i.species_list)+";"+str(i.protein_count)+"_members"]=[i]

# Create a list to guide file writing and sorting it to write each file in a single pass
write_list=list(gbf_dict.keys())
write_list=sorted(write_list,key=lambda x: (x.split(";")[0],int(x.split(";")[1].split("_")[0])))
file_dict={}
for i in write_list:
    try:
        file_dict[i.split(";")[0]].append(i.split(";")[1])
    except KeyError:
        file_dict[i.split(";")[0]]=[i.split(";")[1]]

os.makedirs(args.output_folder+"GenesByFamiliesSpec/",exist_ok=True)
for i in file_dict.keys():
    savefile=open(args.output_folder+"GenesByFamiliesSpec/"+i+".gbfs.csv","w")
    for j in file_dict[i]:
        families=gbf_dict[i+";"+j]
        savefile.write("@ "+j.replace("_"," ")+": "+str(len(families))+" families"+"\n")
        savefile.write(","+i.replace("_",",")+"\n")
        for k in families:
            savefile.write(k.name+",")
            for m in i.split("_"):
                savefile.write(str(len(k.proteins[m]))+",")
            savefile.write("\n")
        savefile.write("\n")
    savefile.close()

# Cleanup temporary lists and dicts
write_list=[]
file_dict={}
gbf_dict={}


print("Creating Fasta files for each groups")
print("This step will take some time and can take a lot of space on computer depending on your groups size")
print("This will also create a lot of folder and files in the Group_Fasta folder of your output")
grps_count=len(ortho_grps)
current_count=0
for group in ortho_grps:
    grp_name=get_valid_filename(group.name)
    
    os.makedirs(args.output_folder+"Groups_Fasta/"+grp_name, exist_ok=True)
    # Create a fasta containing all the proteins of a group
    total_fasta=open(args.output_folder+"Groups_Fasta/"+grp_name+"/"+"allproteins.fasta","w")
    for i in list(group.proteins.items()):
        # Create a species specific fasta
        spc_fasta=open(args.output_folder+"Groups_Fasta/"+grp_name+"/"+i[0]+".fasta","w")
        for j in i[1]:
            total_fasta.write(">"+j+" "+anno_dict[j]+"\n")
            spc_fasta.write(">"+j+" "+anno_dict[j]+"\n")
            total_fasta.write("".join(seq_dict[j])+"\n")
            spc_fasta.write("".join(seq_dict[j])+"\n")
        spc_fasta.close()
    total_fasta.close()
    current_count+=1
    if current_count % 1000 == 0:
        print("%d / %d"%(current_count,grps_count))
print("%d / %d"%(current_count,grps_count))


print("Generation Fasta files according to Ortho Table information")
print("This step will take some time and can take a lot of space on computer depending on your groups size")
print("This will also create a lot of folder and files in the OrthoTable_Fasta folder of your output")
current_count=0
for i in ortho_grps:
    for j in i.species_list:
        if i.classification[0] == "unique":
            save_folders=args.output_folder+"OrthoTable_Fasta/"+j+"/"+i.classification[0]+"_"+j+"/"+i.classification[1]+"/"
        else:
            save_folders=args.output_folder+"OrthoTable_Fasta/"+j+"/"+i.classification[0]+"_"+j+"/"+i.classification[0]+"_"+str(len(i.species_list))+"_species"+"/"+"_".join(i.species_list)+"/"+i.classification[1]+"/"
        os.makedirs(save_folders, exist_ok=True)
        # If singleton, regroup all fasta in one
        if i.protein_count > 1:
            savefile=open(save_folders+get_valid_filename(i.name)+".fasta",'w')
        else:
            savefile=open(save_folders+"singletons.fasta","a")
        for k in list(i.proteins.items()):
            for l in k[1]:
                savefile.write(">"+l+" "+anno_dict[l]+"\n")
                savefile.write("".join(seq_dict[l])+"\n")
        savefile.close()
    current_count+=1
    if current_count % 1000 == 0:
        print("%d / %d"%(current_count,grps_count))
print("%d / %d"%(current_count,grps_count))

# NOW READY FOR ANNOTATION

# Reduce Annotation is optional (if someone don't want to import Spacy for exemple)
class annot_entry:
    def __init__(self,stri,ori_prot_name, thr=args.spacy_thres, mode=0):
        # Mode 0 means using spacy to reduce annotation size
        # Mode 1 means only remove identicals annotations
        self.mode=mode
        if self.mode == 0:
            self.content=nlp(stri)
        else:
            self.content=stri
        self.prot=ori_prot_name # I keep the protein name just in case I need it later
        self.thr=thr
        
    def __eq__(self,othr):
        if self.mode == 0:
            return self.content.similarity(othr.content) >= self.thr
        else:
            return self.content == othr.content
    def __ne__(self,othr):
        if self.mode == 0:
            return self.content.similarity(othr.content) < self.thr
        else:
            return self.content != othr.content
    def __hash__(self):
        # VERY UGLY, BUT IT WORKS...PROBABLY A PERFORMANCE ISSUE...
        # Because the equality is determined by a command during the __eq__ function, it's not possible to use this value as hash
        # If hash are not equal, the equality is not checked
        # That's why I bypass it by assign the same hash value to every entry
        # So it will be possibly be much slower, but I don't see how to do it otherwise
        if self.mode == 0:
            return hash((1))
        else:
            return hash((self.content))
    def __str__(self):
        return str(self.content)
    def __repr__(self):
        return str(self.content)

def regroup_similar(ortho_group, mode=0):
    comp=set()
    original=[]
    for i in list(ortho_group.proteins.values()):
        for j in i:
            comp.add(annot_entry(anno_dict[j],j,mode=mode)) 
    return (ortho_group.name,list(comp)) 

print("Extracting annotation informations from proteomes according to ortholog groups")
if args.ignore_spacy == False:
    print("Spacy will be used to reduce annotation count")
    print("This step will take some time")
results=[]
current_count=0
for i in ortho_grps:
    if i.protein_count==1:
        # If singletons are included in analysis. Only one annotation
        results.append((i.name,[anno_dict[list(i.proteins.values())[0][0]]],list(i.proteins.values())[0][0]))
    else:
        results.append(regroup_similar(i,args.ignore_spacy))
    current_count+=1
    if current_count % 100 == 0:
        print("%d / %d"%(current_count,grps_count))
print("%d / %d"%(current_count,grps_count))    


print("Writing annotation summary TSV")
tsvfile=open(args.output_folder+"annotation_summary.tsv",'w')
tsvfile.write("Group Name\tAnnotation\n")
for i in results:
    # TSV file for each group is not compatible with spacy principle (but I leave some code, just in case)
    for j in i[1]:
        if isinstance(j,annot_entry) :
            tsvfile.write(i[0]+"\t"+str(j)+"\n")
        else:
            tsvfile.write(i[0]+"\t"+str(j)+"\n")

tsvfile.close()

print("\nPROGRAM DONE\n")
