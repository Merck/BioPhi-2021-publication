SHELL := /bin/bash

CONDA_ACTIVATE = source $$(conda info --base)/bin/activate
BIOPHI_ENV_NAME = biophi
CONDA_BIOPHI_ENV = $(CONDA_ACTIVATE) $(BIOPHI_ENV_NAME)

SOURCE_DATA := ../oas-dataset/data
TRAIN_DATA := ../oas-training/data

HUMAN_STUDY_PATHS := $(shell cut -f1 $(SOURCE_DATA)/human/meta/studies.tsv 2>/dev/null | tail -n +2)
MOUSE_STUDY_PATHS := $(shell cut -f1 $(SOURCE_DATA)/mouse/meta/studies.tsv 2>/dev/null | tail -n +2)
ALL_STUDY_PATHS := $(shell cut -f1 $(SOURCE_DATA)/all/meta/studies.tsv 2>/dev/null | tail -n +2)

BIOPHI_OAS_DB := ../../biophi/work/OASis_9mers_v1.db

#-------------#
# Environment #
#-------------#

## Create local Conda environment in "condaenv" folder
condaenv: environment.yml
	hpc/condaenv $<
	hpc/conda-job $@ pip install -e .

#---------------------------#
# Reusable template targets #
#---------------------------#

# Generate heavy chain IMGT numbered CSV using any FASTA file
# Example: Use "make data/my/file_imgt_H.csv" to run ANARCI on "data/my/file_H.fa"
data/%_imgt_H.csv: data/%_H.fa
	anarci -i $< --csv -o $(patsubst %_H.csv,%,$@) --scheme imgt

# Generate light chain IMGT numbered CSV using any FASTA file
# Example: Use "make data/my/file_imgt_KL.csv" to run ANARCI on "make data/my/file_KL.fa"
data/%_imgt_KL.csv: data/%_KL.fa
	anarci -i $< --csv -o $(patsubst %_KL.csv,%,$@) --scheme imgt

# Generate heavy chain Kabat numbered CSV using any FASTA file
# Example: Use "make data/my/file_kabat_H.csv" to run ANARCI on "data/my/file_H.fa"
data/%_kabat_H.csv: data/%_H.fa
	anarci -i $< --csv -o $(patsubst %_H.csv,%,$@) --scheme kabat

# Generate light chain Kabat numbered CSV using any FASTA file
# Example: Use "make data/my/file_kabat_KL.csv" to run ANARCI on "data/my/file_KL.fa"
data/%_kabat_KL.csv: data/%_KL.fa
	anarci -i $< --csv -o $(patsubst %_KL.csv,%,$@) --scheme kabat

# Generate OASis humanness report using any FASTA file
# Example: Use "make data/my/file_oasis.xlsx" to run OASis on "data/my/file.fa"
data/%_oasis.xlsx: data/%.fa
	$(CONDA_BIOPHI_ENV) && biophi oasis \
        $< \
        --output $@ \
        --oasis-db $(BIOPHI_OAS_DB)

# Generate OASis IMGT scheme humanness report using any FASTA file
# Example: Use "make data/my/file_oasis_imgt.xlsx" to run OASis on "data/my/file.fa"
data/%_oasis_imgt.xlsx: data/%.fa
	$(CONDA_BIOPHI_ENV) && biophi oasis \
        $< \
        --output $@ \
        --scheme imgt \
        --oasis-db $(BIOPHI_OAS_DB)

# Generate % germline content using any FASTA file
# Example: Use "make data/my/file_GC.tsv" to run germline content on "make data/my/file.fa"
data/%_GC.tsv: data/%.fa
	bin/humanness_germline_content.py $< $@

# Generate positional germline content using any FASTA file
# Example: Use "make data/my/file_GC_per_position.tsv" to run germline content on "make data/my/file.fa"
data/%_GC_per_position.tsv: data/%.fa
	bin/humanness_germline_content.py $< $@ --scheme kabat --per-position

# Generate T20 humanness score using any FASTA file
# Example: Use "make data/my/file_T20.tsv" to run T20 on "make data/my/file.fa"
data/%_T20.tsv: data/%.fa
	bin/humanness_t20_score.py $< $@

# Generate humanness Z-score using any FASTA file
# Example: Use "make data/my/file_Zscore.tsv" to run Z-score on "data/my/file.fa"
data/%_Zscore.tsv: data/%.fa
	bin/humanness_z_score.py $< $@

# Generate netMHCIIpan predictions using any FASTA file
# Example: Use "make data/my/file_netMHCIIpan.tsv" to run netMHCIIpan on "make data/my/file.fa"
data/%_netMHCIIpan.tsv: data/%.fa
	netMHCIIpan \
        -f $< \
        -length 9 \
        -a DRB1_0101,DRB1_0301,DRB1_0401,DRB1_0701,DRB1_0801,DRB1_1101,DRB1_1301,DRB1_1501 \
        -xls -xlsfile $@

#-----------------#
# Hu-mAb 25 pairs #
#-----------------#

data/tasks/humab_25_pairs/pairs: data/tasks/humab_25_pairs/pairs/experimental_oasis.xlsx data/tasks/humab_25_pairs/pairs/humab_oasis.xlsx data/tasks/humab_25_pairs/pairs/parental_oasis.xlsx data/tasks/humab_25_pairs/pairs/experimental_T20.tsv data/tasks/humab_25_pairs/pairs/humab_T20.tsv data/tasks/humab_25_pairs/pairs/parental_T20.tsv

data/tasks/humab_25_pairs/predicted: data/tasks/humab_25_pairs/predicted/sapiens1_oasis.xlsx data/tasks/humab_25_pairs/predicted/sapiens1_T20.tsv data/tasks/humab_25_pairs/predicted/sapiens2_oasis.xlsx data/tasks/humab_25_pairs/predicted/sapiens2_T20.tsv data/tasks/humab_25_pairs/predicted/sapiens3_oasis.xlsx data/tasks/humab_25_pairs/predicted/sapiens3_T20.tsv data/tasks/humab_25_pairs/predicted/sapiens4_oasis.xlsx data/tasks/humab_25_pairs/predicted/sapiens4_T20.tsv data/tasks/humab_25_pairs/predicted/sapiens5_oasis.xlsx data/tasks/humab_25_pairs/predicted/sapiens5_T20.tsv

data/tasks/humab_25_pairs/predicted/sapiens1.fa: data/tasks/humab_25_pairs/pairs/parental.fa
	$(CONDA_BIOPHI_ENV) && biophi sapiens $< --output $@ --fasta-only --version v1
    
data/tasks/humab_25_pairs/predicted/sapiens2.fa: data/tasks/humab_25_pairs/predicted/sapiens1.fa
	$(CONDA_BIOPHI_ENV) && biophi sapiens $< --output $@ --fasta-only --version v1
    
data/tasks/humab_25_pairs/predicted/sapiens3.fa: data/tasks/humab_25_pairs/predicted/sapiens2.fa
	$(CONDA_BIOPHI_ENV) && biophi sapiens $< --output $@ --fasta-only --version v1

data/tasks/humab_25_pairs/predicted/sapiens4.fa: data/tasks/humab_25_pairs/predicted/sapiens3.fa
	$(CONDA_BIOPHI_ENV) && biophi sapiens $< --output $@ --fasta-only --version v1

data/tasks/humab_25_pairs/predicted/sapiens5.fa: data/tasks/humab_25_pairs/predicted/sapiens4.fa
	$(CONDA_BIOPHI_ENV) && biophi sapiens $< --output $@ --fasta-only --version v1


data/tasks/humab_25_pairs/predicted/sapiens1: data/tasks/humab_25_pairs/pairs/parental.fa
	$(CONDA_BIOPHI_ENV) && biophi sapiens \
        $< \
        --output $@ \
        --version v1 \
        --oasis-db $(BIOPHI_OAS_DB) \
        --iterations 1
        

data/tasks/humab_25_pairs/predicted/sapiens2: data/tasks/humab_25_pairs/pairs/parental.fa
	$(CONDA_BIOPHI_ENV) && biophi sapiens \
        $< \
        --output $@ \
        --version v1 \
        --oasis-db $(BIOPHI_OAS_DB) \
        --iterations 2
        

data/tasks/humab_25_pairs/predicted/sapiens3: data/tasks/humab_25_pairs/pairs/parental.fa
	$(CONDA_BIOPHI_ENV) && biophi sapiens \
        $< \
        --output $@ \
        --version v1 \
        --oasis-db $(BIOPHI_OAS_DB) \
        --iterations 3
        

data/tasks/humab_25_pairs/predicted/sapiens4: data/tasks/humab_25_pairs/pairs/parental.fa
	$(CONDA_BIOPHI_ENV) && biophi sapiens \
        $< \
        --output $@ \
        --version v1 \
        --oasis-db $(BIOPHI_OAS_DB) \
        --iterations 4
        

data/tasks/humab_25_pairs/predicted/sapiens5: data/tasks/humab_25_pairs/pairs/parental.fa
	$(CONDA_BIOPHI_ENV) && biophi sapiens \
        $< \
        --output $@ \
        --version v1 \
        --oasis-db $(BIOPHI_OAS_DB) \
        --iterations 5


data/tasks/humab_25_pairs/predicted/sapiens_chothia1.fa: data/tasks/humab_25_pairs/pairs/parental.fa
	$(CONDA_BIOPHI_ENV) && biophi sapiens $< --output $@ --fasta-only --scheme chothia --cdr-definition chothia --version v1
    
data/tasks/humab_25_pairs/predicted/sapiens_chothia2.fa: data/tasks/humab_25_pairs/predicted/sapiens_chothia1.fa
	$(CONDA_BIOPHI_ENV) && biophi sapiens $< --output $@ --fasta-only --scheme chothia --cdr-definition chothia --version v1
    
data/tasks/humab_25_pairs/predicted/sapiens_chothia3.fa: data/tasks/humab_25_pairs/predicted/sapiens_chothia2.fa
	$(CONDA_BIOPHI_ENV) && biophi sapiens $< --output $@ --fasta-only --scheme chothia --cdr-definition chothia --version v1


#-------------------------------#
# Therapeutic Antibody Metadata #
#-------------------------------#

data/sabdab: data/sabdab/TheraSAbDab.csv

# See http://opig.stats.ox.ac.uk/webapps/newsabdab/therasabdab/search/?all=true
#
## Download therapeutic antibody dataset
data/sabdab/TheraSAbDab_downloaded.csv:
	mkdir -p $(@D)
	wget -O $@ http://opig.stats.ox.ac.uk/webapps/newsabdab/static/downloads/TheraSAbDab_SeqStruc_OnlineDownload.csv

# data/sabdab/TheraSAbDab_clean.tsv: Created in notebook

#-------------------------#
# Therapeutic Rediscovery #
#-------------------------#

data/tasks/therapeutic_rediscovery/oas_hits: data/tasks/therapeutic_rediscovery/oas_hits/heavy data/tasks/therapeutic_rediscovery/oas_hits/light

#
# Look for heavy chains by total identity
#

data/tasks/therapeutic_rediscovery/oas_hits/heavy: $(patsubst %,data/tasks/therapeutic_rediscovery/oas_hits/heavy/%,$(ALL_STUDY_PATHS))

data/tasks/therapeutic_rediscovery/oas_hits/heavy/%: $(SOURCE_DATA)/all/meta/heavy-units-list/%.txt $(SOURCE_DATA)/all/json
	mkdir -p $@
	@if [ -s $< ]; then \
				printf "make -j 32" > $@/run.sh; \
				while read name; do \
						printf " $@/$$name.csv" >> $@/run.sh; \
				done < $<; \
				hpc/conda-job $@ bash $@/run.sh; \
	else \
				echo "Creating empty dir for $*, no units in $<"; \
	fi

data/tasks/therapeutic_rediscovery/oas_hits/heavy/%.csv: $(SOURCE_DATA)/all/json/%.json.gz data/tasks/therapeutic_rediscovery/thera/humanized_imgt_H.csv
	mkdir -p $(@D)
	bin/global_search_imgt_oas.py \
				$< \
				--query $(word 2,$^) \
				--output $@

#
# Look for light chains by total identity
#

data/tasks/therapeutic_rediscovery/oas_hits/light: $(patsubst %,data/tasks/therapeutic_rediscovery/oas_hits/light/%,$(ALL_STUDY_PATHS))

data/tasks/therapeutic_rediscovery/oas_hits/light/%: $(SOURCE_DATA)/all/meta/light-units-list/%.txt $(SOURCE_DATA)/all/json
	mkdir -p $@
	@if [ -s $< ]; then \
				printf "make -j 32" > $@/run.sh; \
				while read name; do \
						printf " $@/$$name.csv" >> $@/run.sh; \
				done < $<; \
				hpc/conda-job $@ bash $@/run.sh; \
	else \
				echo "Creating empty dir for $*, no units in $<"; \
	fi

data/tasks/therapeutic_rediscovery/oas_hits/light/%.csv: $(SOURCE_DATA)/all/json/%.json.gz data/tasks/therapeutic_rediscovery/thera/humanized_imgt_KL.csv
	mkdir -p $(@D)
	bin/global_search_imgt_oas.py \
				$< \
				--query $(word 2,$^) \
				--output $@




data/tasks/therapeutic_rediscovery/oas_cdr_hits: data/tasks/therapeutic_rediscovery/oas_cdr_hits/heavy data/tasks/therapeutic_rediscovery/oas_cdr_hits/light

#
# Look for heavy chains by prioritizing CDR identity
#

data/tasks/therapeutic_rediscovery/oas_cdr_hits/heavy: $(patsubst %,data/tasks/therapeutic_rediscovery/oas_cdr_hits/heavy/%,$(ALL_STUDY_PATHS))

data/tasks/therapeutic_rediscovery/oas_cdr_hits/heavy/%: $(SOURCE_DATA)/all/meta/heavy-units-list/%.txt $(SOURCE_DATA)/all/json
	mkdir -p $@
	@if [ -s $< ]; then \
				printf "make -j 16" > $@/run.sh; \
				while read name; do \
						printf " $@/$$name.csv" >> $@/run.sh; \
				done < $<; \
				hpc/conda-job $@ bash $@/run.sh; \
	else \
				echo "Creating empty dir for $*, no units in $<"; \
	fi

data/tasks/therapeutic_rediscovery/oas_cdr_hits/heavy/%.csv: $(SOURCE_DATA)/all/json/%.json.gz data/tasks/therapeutic_rediscovery/thera/humanized_imgt_H.csv
	mkdir -p $(@D)
	bin/cdr_search_imgt_oas.py \
				$< \
				--query $(word 2,$^) \
				--output $@

#
# Look for light chains by prioritizing CDR identity
#

data/tasks/therapeutic_rediscovery/oas_cdr_hits/light: $(patsubst %,data/tasks/therapeutic_rediscovery/oas_cdr_hits/light/%,$(ALL_STUDY_PATHS))

data/tasks/therapeutic_rediscovery/oas_cdr_hits/light/%: $(SOURCE_DATA)/all/meta/light-units-list/%.txt $(SOURCE_DATA)/all/json
	mkdir -p $@
	@if [ -s $< ]; then \
				printf "make -j 16" > $@/run.sh; \
				while read name; do \
						printf " $@/$$name.csv" >> $@/run.sh; \
				done < $<; \
				hpc/conda-job $@ bash $@/run.sh; \
	else \
				echo "Creating empty dir for $*, no units in $<"; \
	fi

data/tasks/therapeutic_rediscovery/oas_cdr_hits/light/%.csv: $(SOURCE_DATA)/all/json/%.json.gz data/tasks/therapeutic_rediscovery/thera/humanized_imgt_KL.csv
	mkdir -p $(@D)
	bin/cdr_search_imgt_oas.py \
				$< \
				--query $(word 2,$^) \
				--output $@



#
# Generate BioPhi Sapiens humanized sequences
#

# Automatic germline
    
data/tasks/therapeutic_rediscovery/predicted_automatic: data/tasks/therapeutic_rediscovery/predicted_automatic/sapiens1_oasis.xlsx data/tasks/therapeutic_rediscovery/predicted_automatic/sapiens1_T20.tsv

data/tasks/therapeutic_rediscovery/predicted_automatic/sapiens1.fa: data/tasks/therapeutic_rediscovery/oas_cdr_hits/parental.fa
	$(CONDA_BIOPHI_ENV) && biophi sapiens $< --output $@ --fasta-only --version v1
    
# Manual germline    

data/tasks/therapeutic_rediscovery/predicted_manual: data/tasks/therapeutic_rediscovery/predicted_manual/sapiens1_oasis.xlsx data/tasks/therapeutic_rediscovery/predicted_manual/sapiens1_T20.tsv

data/tasks/therapeutic_rediscovery/predicted_manual/sapiens1.fa: data/tasks/therapeutic_rediscovery/cdr_grafts_manual/vernier_grafts.fa
	$(CONDA_BIOPHI_ENV) && biophi sapiens $< --output $@ --fasta-only --version v1
    
    

#-------------#
# Deamidation #
#-------------#

data/tasks/ptm_mitigation/sapiens: data/tasks/ptm_mitigation/ptm_seqs_masked.fa
	$(CONDA_BIOPHI_ENV) && biophi sapiens \
        $< \
        --output $@ \
        --version v1 \
        --oasis-db $(BIOPHI_OAS_DB) \
        --iterations 1

