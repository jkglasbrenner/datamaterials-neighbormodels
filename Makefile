SHELL				=	/bin/sh

PROJNAME			=	neighbormodels

RM					=	rm
ECHO				=	echo
COPY				=	cp
FIND				=	find

CONDA				=	conda
CONDA_ENV_FILE		=	environment.yaml

PYTHON				=	/usr/bin/env python
PYTHON_SETUP		=	setup.py
PYTHON_SETUP_DOCS	=	build_sphinx

ALL_FILES			=

CLEAN_FILES			=	*_cache/												\
						docs/_build/*

define cleanup
	-$(RM) -rf $(CLEAN_FILES)
endef

define pycache_cleanup
	$(FIND) -name "__pycache__" -type d -exec $(RM) -rf {} +
endef

define setup_environment
	bash -lc "$(CONDA) env update --file $(CONDA_ENV_FILE)"
endef

define make_docs
	$(PYTHON) ./$(PYTHON_SETUP) $(PYTHON_SETUP_DOCS)
endef

.SILENT		:
.PHONY		:	all clean docs environment

all			:	$(ALL_FILES)

docs		:
	$(call make_docs)

clean		:
	$(call cleanup)
	$(call pycache_cleanup)

environment	:
	$(ECHO) Setting up conda environment
	$(call setup_environment)
