SHELL				=	/bin/sh

RM					=	rm
ECHO				=	echo
COPY				=	cp

CONDA				=	conda
CONDA_ENV_FILE		=	environment.yaml
CONDA_ENV_NAME		=	crystal-structures
CONDA_ENV_ACTIVATE	=	$(CONDA) activate

ALL_FILES			=

CLEAN_FILES	=	*_cache/

define cleanup
	-$(RM) -rf $(CLEAN_FILES)
endef

define setup_environment
	bash -lc "$(CONDA) env update --file $(CONDA_ENV_FILE)"
endef

.SILENT	:
.PHONY	:	all clean environment

all	: $(ALL_FILES)

clean	:
	$(call cleanup)

environment	:
	$(ECHO) Setting up conda environment
	$(call setup_environment)
