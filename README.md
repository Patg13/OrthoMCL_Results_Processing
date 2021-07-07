# OrthoMCL_Results_Processing (OMRP)
This Python3 program will take the orthomcl output (groups, proteomes and singletons) and produce different human readable results

# Python3 Dependencies
**pyvenn** (_already available in this repo as venn_l_, original reference: https://github.com/tctianchi/pyvenn)

**matplotlib 3.3.4 (Required)** -> Used by pyvenn to generate venn diagrams

**Spacy 3.0.3 (Optional, use --ignore_spacy flag to skip)** -> Regroup similar protein function annotation according to their textual similarity

**en_core_web_lg dictionnary for Spacy** -> Gives best results during testing

# Usage:
  
  -g GROUP_FILE, --group GROUP_FILE
  
                        Group file containing proteins groups [REQUIRED]
                        
  -s SINGLETONS_FILE, --singleton SINGLETONS_FILE
  
                        Singleton file containing single proteins [OPTIONAL]
                        
  -p PROTEOME_FILE, --proteome PROTEOME_FILE
  
                        Fasta file regrouping all proteomes (MUST INCLUDE
                        PROTEIN ANNOTATION AND ID MUST BE SEPERATE BY A SPACE)
                        [OPTIONAL]
                        
  -c SPC_LIST, --species SPC_LIST
  
                        Species tag list used in data (ex: SP1,SP2,...) [REQUIRED]
                        
  -r REGROUP_LIST, --regroup REGROUP_LIST
  
                        List of species you want to regroup (; to link ,to
                        seperate) SP1;SP2,SP3;SP4 means 1 regroup with SP1 and
                        SP2 and 1 regroup with SP3 and SP4 [OPTIONAL]
                        
  -l TABLE_LIMIT, --limit TABLE_LIMIT
  
                        Ortho Table limit for cluster count (every cluster >=
                        this value will have their count add to the same value
                        (Val+)) ! higher value means more columns in the table
                        ! [DEFAULT=11]
                        
  --venn_dpi DPI
  
                        Pixel Density for Venn Diagram greater value means
                        sharper but bigger png [DEFAULT=400]
                        
  --ignore_spacy
  
                        Ignore Spacy step to reduce annotation count (remove
                        spacy dependency, MUCH faster, but will produce a
                        massive annotation CSV)
                        
  --spacy_simthres SPACY_THRES
  
                        Similarity threshold for spacy annotation regrouping
                        MODIFY AT YOUR OWN RISK [DEFAULT=0.95]
                        
  --spacy_pipeline SPACY_PIPE
  
                        Spacy language pipeline to use (all the developpement was
                        done with en_core_web_lg, so it's not recommended to
                        use another) [default=en_core_web_lg]
                        
  -o OUTPUT_FOLDER, --output OUTPUT_FOLDER
  
                        Output folder name [Default=orthomcl_resgen_out/]
                        
                        
 # Input file
 
This program takes the OrthoMCL output files (singletons.txt and groups.txt) as is. The proteomes file must be generated by the user (a simple cat *.fasta > proteomes.fasta will do the trick) DO NOT USE THE COMPLIANT FASTA SEQUENCE GENERATED BY ORTHOMCL, it must be the original proteome files
 
All the proteomes sequences must have a unique sequence identifier AND a protein functional annotation
 
The singletons file can be ommited if you don't want to include the singletons in your results

# OrthoTable
After preparing the groups according to the species and regroups lists, OMRP will compile the protein counts into the OrthoTable (which contain the group count for every species combination possible andn their protein count (up to --limit). This table will be save in OrthoTable.tsv

# Venn Diagrams
Using the information in OrthoTable, 3 venn diagrams will be generated; B C and P. P_venn will use the protein count, C will use the group count and B will contain both information (group_count on top and protein count in parenthesis at the bottom). Please note that the program only support up to 6 species for the venn diagram, if you have more then that, you must use the --regroup option the regroup some species together (or you can leave it like that, but venn diagram generation will be skipped). Also note that a level 6 venn diagram is very hard to read (using triangles instead of ellipses).

# Spacy and annotation extraction
Using the proteome file, OMRP will create a TSV file containing all annotation for each group identified by OrthoMCL. Because there can be a lot of different annotation (depending on the database you use), OMRP uses the SPACY library to detect similar annotation according to their textual similarity (exemple, DNase Acc 3344 and DNase Acc 6677 will be considered the same annotation and only one will appear in the TSV). You can completely ignore the spacy step, which will result in only identical annotation will be regrouped (the TSV will be much bigger though). 

# Genes By Families by Species
This step will generate csv files for each shared line in the OrthoTable according to this syntax (in this exemple, the csv filename would be **SP1_SP2.gpfs.csv**):

@ 3 members: 23 families

	SP1	SP2
	
grp_8045:	1	2

grp_8046:	1	2

@ 4 members: 6 families	

	SP1	SP2
	
grp_4165:	1	3


This is very useful to determine the species gene count for each Ortho Group 

# Groups Fasta files
OMRP will generate a serie of folder (one for each ortho group) which will contain multiple fasta files (one file containing all proteins for this group and one file for each species). Depending on your ortho group count, this can take a lot of space and will generate a lot of files. If you included a singletons file, they will be noted as single_# instead of your chosen prefix.

# OrthoTable Fasta
ORMP will also generate another serie of folder which will reflect the content of the OrthoTable. Each fasta will be classified by their respective base species, then if they are in a shared state (ex : shared between species 1 and species 2) or a unique state (only in the current species), then by the specific share (specie1_species2), then by the number of protein shared (2 - TABLE_LIMIT+ ), uou will then find every Ortho group fasta. This is probably the most heavy step for space and file count.
